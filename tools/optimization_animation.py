"""Pack accepted ODFR states and render them without loading the history into RAM.

Visual specification: 960 x 600 white canvas, original experiment (gray dashed)
and accepted model (orange solid), shared fixed axes over the entire history,
no interpolation/smoothing. SSE_W is the original full-fit objective, not an
uncertainty estimate or a loss recomputed from rounded visualization values.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import warnings
import zipfile

import numpy as np


def read_table(path: Path):
    text = path.read_text(encoding="ascii")
    if not text.rstrip().endswith("# complete"):
        raise ValueError(f"Incomplete state table: {path}")
    metadata = {}
    for line in text.splitlines():
        if line.startswith("# ") and line != "# complete":
            key, value = line[2:].split(maxsplit=1)
            if key != "spectrum":
                metadata[key] = value
    data = np.loadtxt(io.StringIO(text), ndmin=2)
    if data.shape[1] != 3 or not np.isfinite(data).all():
        raise ValueError(f"Invalid state table: {path}")
    return metadata, data


class History:
    """Reiterable state stream: one full table plus one delta at a time."""

    def __init__(self, source):
        self.source = Path(source)
        self.relative_precision = None
        if self.source.is_dir():
            self.files = sorted(self.source.glob("step_*.dat"))
            if not self.files:
                raise ValueError("No step_*.dat files found")
            # A still-running fit may have one unfinished trailing file.
            try:
                read_table(self.files[-1])
            except ValueError:
                warnings.warn(f"Ignoring incomplete trailing state: {self.files[-1]}")
                self.files.pop()
            if not self.files:
                raise ValueError("No complete state tables found")
            self.experiment = np.loadtxt(self.source / "experiment.dat", ndmin=2)
            self.count = len(self.files)
        else:
            with zipfile.ZipFile(self.source) as archive:
                self.manifest = json.loads(archive.read("manifest.json"))
                if self.manifest["format"] != "odfr-steps-v1":
                    raise ValueError("Unsupported archive format")
                self.experiment = np.load(io.BytesIO(archive.read("experiment.npy")), allow_pickle=False)
            self.relative_precision = self.manifest["relative_precision"]
            self.count = len(self.manifest["states"])

    def __iter__(self):
        if self.source.is_dir():
            for path in self.files:
                yield read_table(path)
            return
        with zipfile.ZipFile(self.source) as archive:
            def load(name):
                return np.load(io.BytesIO(archive.read(name)), allow_pickle=False)
            grid = load("grid.npy")
            base = load("frames/000000.npy")
            offsets = np.zeros(len(base), dtype=np.int64)
            quantum = load("quantum.npy") if self.relative_precision else None
            for i, meta in enumerate(self.manifest["states"]):
                payload = base if i == 0 else load(f"frames/{i:06d}.npy")
                if i and quantum is not None:
                    offsets += payload
                    model = base + offsets * quantum
                else:
                    model = payload
                yield meta, np.column_stack((grid, model))


def write_array(archive, name, array):
    with archive.open(name, "w") as member:
        np.lib.format.write_array(member, np.asarray(array), allow_pickle=False)


def pack(source, output, relative_precision=None):
    history = History(source)
    if not history.source.is_dir():
        raise ValueError("Pack from original tables, not an already rounded archive")
    if relative_precision is not None and not 1e-12 <= relative_precision <= 1:
        raise ValueError("--relative-precision must be between 1e-12 and 1")
    _, initial = next(iter(history))
    grid, base = initial[:, :2].copy(), initial[:, 2].copy()
    ids = sorted({int(i) for i in grid[:, 0]})
    peaks = {i: 0.0 for i in ids}
    for _, data in history:
        if not np.array_equal(data[:, :2], grid):
            raise ValueError("The spectral grids changed; keep the original tables")
        for i in ids:
            peaks[i] = max(peaks[i], float(np.max(np.abs(data[data[:, 0] == i, 2]))))
    quantum = np.array([(peaks[int(i)] or 1.0) * (relative_precision or 0) for i in grid[:, 0]])
    states, max_errors = [], {i: 0.0 for i in ids}
    previous = np.zeros(len(base), dtype=np.int64)
    # Exclusive creation: never replace an archive or modify its source tables.
    with zipfile.ZipFile(output, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        write_array(archive, "grid.npy", grid)
        write_array(archive, "experiment.npy", history.experiment)
        if relative_precision:
            write_array(archive, "quantum.npy", quantum)
        for index, (meta, data) in enumerate(history):
            model = data[:, 2]
            payload = model.astype("<f8")
            if index and relative_precision:
                # Quantize offsets from the exact first curve, NOT successive
                # floating differences. Integer telescoping prevents drift.
                offsets = np.rint((model - base) / quantum).astype(np.int64)
                restored = base + offsets * quantum
                error = np.abs(restored - model)
                bound = quantum / 2 + 8 * np.finfo(float).eps * np.maximum(np.abs(base), np.abs(model))
                if np.any(error > bound):
                    raise ValueError("Quantization error bound was exceeded")
                for i in ids:
                    max_errors[i] = max(max_errors[i], float(np.max(error[grid[:, 0] == i])))
                payload = (offsets - previous).astype("<i8")
                previous = offsets
            write_array(archive, f"frames/{index:06d}.npy", payload)
            states.append(meta)
        manifest = dict(format="odfr-steps-v1", relative_precision=relative_precision,
                        encoding="integer-offset-deltas" if relative_precision else "float64",
                        spectrum_peak_amplitude=peaks, measured_max_absolute_error=max_errors,
                        states=states, source=str(history.source.resolve()),
                        loss_semantics="Original exact full-fit SSE; never recomputed from rounded curves")
        archive.writestr("manifest.json", json.dumps(manifest, indent=2))
    print(f"Saved {len(states)} states: {output} ({Path(output).stat().st_size:,} bytes)")
    if relative_precision:
        print("Max measured absolute error per spectrum:", max_errors)


def render(source, output, spectrum=1, fps=8, field_unit="ед. входных данных"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from PIL import GifImagePlugin, Image

    if not 0 < fps <= 100:
        raise ValueError("--fps must be in (0, 100]")
    initial_preview = Path(output).with_name(Path(output).stem + "_initial.png")
    if Path(output).suffix.lower() != ".gif":
        raise ValueError("The output filename must end in .gif")
    history = History(source)
    exp = history.experiment[history.experiment[:, 0] == spectrum]
    if not len(exp):
        raise ValueError(f"Spectrum {spectrum} is absent")
    xmin, xmax = np.min(exp[:, 1]), np.max(exp[:, 1])
    ymin, ymax = np.min(exp[:, 2]), np.max(exp[:, 2])
    for _, data in history:
        curve = data[data[:, 0] == spectrum]
        if not len(curve):
            raise ValueError(f"Spectrum {spectrum} is absent from a state")
        xmin, xmax = min(xmin, np.min(curve[:, 1])), max(xmax, np.max(curve[:, 1]))
        ymin, ymax = min(ymin, np.min(curve[:, 2])), max(ymax, np.max(curve[:, 2]))
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 12,
                         "text.color": "#171717", "axes.labelcolor": "#171717",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.7, "svg.fonttype": "none"}):
        fig, ax = plt.subplots(figsize=(9.6, 6), dpi=100)
        fig.subplots_adjust(left=.12, right=.97, bottom=.14, top=.73)
        fig.text(.12, .94, f"Подгонка спектра {spectrum}", fontsize=18, weight="semibold")
        subtitle = fig.text(.12, .875, "", fontsize=12)
        if history.relative_precision:
            fig.text(.12, .825, f"Округление для показа: ≤ {history.relative_precision / 2:g} от пика каждого спектра",
                     color="#5E6462", fontsize=10)
        ax.plot(exp[:, 1], exp[:, 2], color="#747976", lw=1.2, ls="--", label="Эксперимент")
        line, = ax.plot([], [], color="#E67E22", lw=1.8, label="Модель")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, frameon=False, borderaxespad=0)
        padx, pady = max((xmax-xmin)*.02, 1e-12), max((ymax-ymin)*.08, 1e-12)
        ax.set(xlim=(3320, 3410), ylim=(-0.5, 0.6),
               xlabel=f"Магнитное поле ({field_unit})", ylabel="Интенсивность (отн. ед.)")
        ax.xaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.grid(axis="both", color="#D5D8D6", linewidth=.6)
        ax.set_axisbelow(True)
        # Pillow's frame encoder writes one palette/frame at a time. Its normal
        # save_all/PillowWriter path can retain every raster in memory.
        with Path(output).open("wb") as stream:
            for index, (meta, data) in enumerate(history):
                curve = data[data[:, 0] == spectrum]
                line.set_data(curve[:, 1], curve[:, 2])
                stage = "Начальное состояние" if index == 0 else f"Принятое состояние {index}"
                subtitle.set_text(f"{stage} · итерация {meta['iteration']} · SSE_W = {float(meta['SSE_W']):.4g}")
                if index == 0:
                    with initial_preview.open("wb") as preview_stream:
                        fig.savefig(preview_stream, format="png", dpi=200)
                fig.canvas.draw()
                raster = Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert("RGB")
                palette = raster.quantize(colors=256)
                if index == 0:
                    for block in GifImagePlugin.getheader(palette, info={"loop": 0})[0]:
                        stream.write(block)
                duration = max(10, round(100 / fps) * 10)
                for block in GifImagePlugin.getdata(palette, duration=duration, disposal=2,
                                                    include_color_table=True):
                    stream.write(block)
            stream.write(b";")
        # A vector final frame and PNG permit inspection at the delivery size.
        for suffix in (".png", ".svg"):
            preview = Path(output).with_suffix(suffix)
            with preview.open("wb") as stream:
                fig.savefig(stream, format=suffix[1:], dpi=200 if suffix == ".png" else 100)
        plt.close(fig)
    print(f"Rendered {history.count} states: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    compact = commands.add_parser("pack", help="Create a separate compressed archive; preserve all tables")
    compact.add_argument("source", type=Path)
    compact.add_argument("output", type=Path)
    compact.add_argument("--relative-precision", type=float, help="Optional lossy quantum / per-spectrum peak (e.g. 1e-6)")
    draw = commands.add_parser("render", help="Render tables or archive to GIF, initial PNG, final PNG and SVG")
    draw.add_argument("source", type=Path)
    draw.add_argument("output", type=Path)
    draw.add_argument("--spectrum", type=int, default=1)
    draw.add_argument("--fps", type=float, default=8)
    draw.add_argument("--field-unit", default="ед. входных данных")
    args = vars(parser.parse_args())
    command = args.pop("command")
    (pack if command == "pack" else render)(**args)


if __name__ == "__main__":
    main()
