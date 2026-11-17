"""Exercise ODFR through namODFR.nam with isolated copies of smoke inputs."""

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def prepare(case_dir, inputs, jobs, output="out.odf", *, loose_cg_jobs=(), orientations=40):
    """Each job is a list of spectrum variants and produces independent outputs."""
    case_dir.mkdir()
    shutil.copyfile(inputs / "smoke.esr", case_dir / "exp.esr")
    odf = (inputs / "ODFR4F.par").read_text(encoding="ascii").splitlines()
    odf[24] = "20 0"  # Request graphical coefficients for the second site, too.
    (case_dir / "odf.par").write_text("\n".join(odf) + "\n", encoding="ascii")
    template = (inputs / "ODFR4.par").read_text(encoding="ascii").splitlines()
    rows = []
    for job, variants in enumerate(jobs, 1):
        odf_output = output if job == 1 else f"out{job}.odf"
        quoted = "'" + odf_output.replace("'", "''") + "'"
        rows.extend([str(job), f"odf.par,{quoted},{len(variants)}"])
        for index, variant in enumerate(variants):
            name = f"{job}-{index}"
            lines = template.copy()
            lines[0] = "2,0,2"  # Single calculation.
            lines[4] = "57.0,1.0,0"  # Unit weights.
            lines[31] = "0"  # Two independent sites.
            for row in (9, 34):
                lines[row] = f"{orientations},0"  # NORT for each site.
            for row in (13, 38):
                lines[row] = "4,3,2,2,2"  # Small basis for this test.
            lines[17] = f"{1.5 + variant * 0.07},10,0.6,0,0,0,0,0"
            for row in (30, 55):
                lines[row] = f"0,0,{variant * 7},0,0,0"
            lines[19] = f"{2.0087 + variant * 0.001},10,2.006,0,2.0023,0"
            if job in loose_cg_jobs:
                # Later spectra inherit CG tolerance from the first PAR of the job.
                lines[12], lines[37] = "0.15,1000", "0.15,500"
            (case_dir / f"p{name}.par").write_text(
                "\n".join(lines) + "\n", encoding="ascii"
            )
            rows.append(f"p{name}.par,r{name}.par,exp.esr,o{name}.dat,2")
    (case_dir / "namODFR.nam").write_text(
        "\n".join(rows + ["end"]) + "\n", encoding="ascii"
    )
    return case_dir


def run(case_dir, command, *, jobs=1, error_path=None):
    proc = subprocess.run(
        command, cwd=case_dir, capture_output=True, text=True, errors="replace",
        timeout=180, env={**os.environ, "ODFR_OPTIMIZATION_STEPS": "0"},
    )
    if error_path is None:
        if proc.returncode == 0 and proc.stdout.count("Final loss reporting:") == jobs:
            return
    elif proc.returncode and f"Cannot write ODF output {error_path}:" in proc.stdout + proc.stderr:
        return
    raise AssertionError(
        f"{case_dir.name}: exit={proc.returncode}, expected error path={error_path!r}\n"
        f"{proc.stdout}\n{proc.stderr}"
    )


def spectrum(path):
    result = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if len(fields) == 6:
            try:
                # Model field, combined spectrum and contributions of both sites.
                result.append(tuple(float(value.replace("D", "E")) for value in fields[2:]))
            except ValueError:
                continue
    if len(result) != 512 or any(not math.isfinite(x) for row in result for x in row):
        raise AssertionError(f"{path.name}: expected 512 finite model rows")
    return result


def compare(path, reference):
    for index, (actual, expected) in enumerate(zip(spectrum(path), reference), 1):
        # The DAT format rounds model values before writing.
        if any(not math.isclose(a, b, rel_tol=3e-5, abs_tol=1e-10)
               for a, b in zip(actual, expected)):
            raise AssertionError(f"{path.name}, model row {index}: {actual} != {expected}")


def check_batch(root, inputs, command):
    variants = [0, 1, 2] * 4 + [0, 1]
    # 14 spectra * 2 independent sites * 40 orientations exceeds the old 800
    # slots. A second job in the same process checks state reset.
    jobs = [variants, [2, 1]]
    compare_jobs(root, inputs, command, jobs, loose_cg_jobs=(2,))


def check_series(root, inputs, command):
    # MaxCALC permits ten jobs per process. Alternate spectrum count and CG
    # tolerance with a small orientation grid to keep repeated runs affordable.
    compare_jobs(root, inputs, command, [[0, 1], [2]] * 5,
                 loose_cg_jobs=(2, 4, 6, 8, 10), orientations=4)


def compare_jobs(root, inputs, command, jobs, *, loose_cg_jobs, orientations=40):
    batch = prepare(root / "batch", inputs, jobs, loose_cg_jobs=loose_cg_jobs,
                    orientations=orientations)
    run(batch, command, jobs=len(jobs))
    references = {}
    for job, sequence in enumerate(jobs, 1):
        loose = job in loose_cg_jobs
        for variant in sorted(set(sequence)):
            key = (variant, loose)
            if key not in references:
                single = prepare(
                    root / f"single-{variant}-{int(loose)}", inputs, [[variant]],
                    loose_cg_jobs=(1,) if loose else (), orientations=orientations,
                )
                run(single, command)
                references[key] = spectrum(single / "o1-0.dat")
        refs = [references[variant, loose] for variant in sorted(set(sequence))]
        if len({tuple(rows) for rows in refs}) != len(refs):
            raise AssertionError("Expected distinct isolated reference spectra")
        for index, variant in enumerate(sequence):
            compare(batch / f"o{job}-{index}.dat", references[variant, loose])


def output_bytes(case_dir, name):
    paths = [name, "G1" + name, "G2" + name, "o1-0.dat"]
    outputs = [(case_dir / path).read_bytes() for path in paths]
    if not all(outputs):
        raise AssertionError(f"{case_dir.name}: one of {paths} is empty")
    return outputs


def check_output_paths(root, inputs, command):
    flat = prepare(root / "flat", inputs, [[0]])
    run(flat, command)
    expected = output_bytes(flat, "out.odf")

    output = "ODF files/sub/out.odf"
    nested = prepare(root / "nested", inputs, [[0]], output)
    for _ in range(2):  # Missing directories first, existing output on rerun.
        run(nested, command)
        if output_bytes(nested, output) != expected:
            raise AssertionError("Creating parent directories changed the outputs")

    blocked = prepare(root / "blocked-parent", inputs, [[0]], "blocked/out.odf")
    marker = blocked / "G1blocked"
    marker.write_bytes(b"preserve")
    run(blocked, command, error_path="G1blocked/out.odf")
    if marker.read_bytes() != b"preserve":
        raise AssertionError("The file blocking a parent directory was modified")

    blocked = prepare(root / "blocked-open", inputs, [[0]])
    target = blocked / "G2out.odf"
    target.mkdir()
    run(blocked, command, error_path=target.name)
    if not target.is_dir():
        raise AssertionError("The directory blocking an output file was replaced")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    checks = {"batch": check_batch, "output_paths": check_output_paths, "series": check_series}
    parser.add_argument("scenario", choices=checks)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.command:
        parser.error("expected an executable or MPI launch command")
    args.run_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=args.scenario + "-", dir=args.run_root) as temp:
        checks[args.scenario](Path(temp), args.inputs, args.command)


if __name__ == "__main__":
    main()
