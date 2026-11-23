"""Exercise the cross-evaluation cache during fitting and across jobs."""

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


ORIENTATIONS = 3
NPOINTS = 65
STAT_LINE = re.compile(
    r"ORIENTATION_CACHE rank=(\d+) full=(\d+) cfonly=(\d+) spectra=(\d+)"
    r" cross=(\d+) eval=(\d+)"
)


def fit_params(source):
    params = (source / "ODFR4.par").read_text(encoding="ascii").splitlines()[:31]
    params[0] = "3 0 2"  # NL2SOL, with covariance.
    params[1] = "2 1"  # Two varied EPR parameters, one site.
    params[2] = f"{NPOINTS} 1"
    params[4] = "0.0 1.0 0"  # Unit weights.
    params[5] = "3300 -64 64"  # An integer field grid for the synthetic data.
    params[6] = "1e-3 1e-3 1e-5"
    params[9] = f"{ORIENTATIONS} 0"
    params[12] = "0.001 500"
    params[13] = "6 5 4 3 2"
    params[17] = "1.5 1  0.6 0  0.0 0  0.0 0"  # Gaussian, CFS/convolution only.
    params[22] = "4.0E+06 2  4.0E+06 4  4.0E+08 0"  # Diffusion, new CG matrix.
    params[30] = "0.0 0  30.0 0  0.0 0"
    return params


def write_odf(root, coefficient):
    (root / "odf.par").write_text(
        f"0\n0 0.0 80\nODF\n2 0\nCoefficients\n1.0 0\n{coefficient} 1 0.0 0\n",
        encoding="ascii",
    )


def make_experiment(root, source, command):
    root.mkdir()
    shutil.copyfile(source / "smoke.esr", root / "smoke.esr")
    params = fit_params(source)
    params[0] = "2 0 2"
    params[17] = "1.6 0  0.6 0  0.0 0  0.0 0"
    params[22] = "4.3E+06 0  4.3E+06 0  4.0E+08 0"
    (root / "seed.par").write_text("\n".join(params) + "\n", encoding="ascii")
    write_odf(root, 0.32)
    (root / "namODFR.nam").write_text(
        "1\nodf.par, seed-out.par, 1\nseed.par, r0.par, smoke.esr, r0.dat, 2\nend\n",
        encoding="ascii",
    )
    run(root, command, "0")
    lines = (root / "r0.dat").read_text(encoding="ascii").splitlines()
    first = next(i + 1 for i, line in enumerate(lines) if "Hmod" in line)
    # Keep the model columns and extend the grid by one point on either end
    # so the interpolation in choiceExp covers the full fitting interval.
    samples = [line.split()[2:4] for line in lines[first:first + NPOINTS]]
    assert len(samples) == NPOINTS and all(len(sample) == 2 for sample in samples), samples
    result = root / "fit.esr"
    result.write_text("Spectr:\n3234 0\n" +
                      "".join(f"{field} {model}\n" for field, model in samples) +
                      "3366 0\n", encoding="ascii")
    return result


def prepare(root, source, experiment):
    root.mkdir()
    shutil.copyfile(experiment, root / "fit.esr")
    params = fit_params(source)
    (root / "nl2.par").write_text("\n".join(params) + "\n", encoding="ascii")
    params[0] = "4 0 2"  # MINPACK, including its legacy Jacobian path.
    (root / "minpack.par").write_text("\n".join(params) + "\n", encoding="ascii")
    write_odf(root, 0.3)
    rows = []
    for task, engine in enumerate(("nl2", "minpack", "nl2"), 1):
        rows.extend(
            (str(task), f"odf.par, o{task}.par, 1",
             f"{engine}.par, r{task}.par, fit.esr, r{task}.dat, 2")
        )
    (root / "namODFR.nam").write_text("\n".join(rows + ["end"]) + "\n", encoding="ascii")


def run(root, command, mode):
    env = dict(os.environ, ODFR4_ORIENTATION_CACHE=mode, ODFR_OPTIMIZATION_STEPS="1")
    result = subprocess.run(command, cwd=root, env=env, timeout=90, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert result.returncode == 0, result.stdout
    return result.stdout


def summarize(stats, ranks, mode):
    """Sum ranks at each evaluation without relying on MPI stdout ordering."""
    per_rank = {rank: [] for rank in range(ranks)}
    for rank, full, cfonly, spectra, cross, evaluation in STAT_LINE.findall(stats):
        jobs = per_rank[int(rank)]
        evaluation = int(evaluation)
        if evaluation == 1:
            jobs.append([])
        assert jobs and evaluation == len(jobs[-1]) + 1, (mode, rank, evaluation, jobs)
        jobs[-1].append(tuple(map(int, (full, cfonly, spectra, cross))))

    # The two NL2 jobs report every evaluation. MINPACK runs without cache stats.
    for rank, jobs in per_rank.items():
        assert len(jobs) == 2 and all(len(job) >= 4 for job in jobs), (mode, rank, jobs)
    lengths = [[len(job) for job in jobs] for jobs in per_rank.values()]
    assert all(length == lengths[0] for length in lengths), (mode, lengths)
    return [
        [tuple(map(sum, zip(*evaluation))) for evaluation in zip(*job)]
        for job in zip(*per_rank.values())
    ]


def read_outputs(root):
    names = [f"{kind}{task}.{ext}" for task in range(1, 4)
             for kind, ext in (("o", "par"), ("G1o", "par"),
                               ("r", "par"), ("r", "dat"))]
    names += ["temp 1.dat", "temp 1.par", "tempODF.par", "G1tempODF.par"]
    outputs = {name: (root / name).read_bytes() for name in names}
    assert all(outputs.values()), names
    fits = sorted((root / "optimization_steps").glob("fit_*"))
    assert len(fits) == 2, fits
    for fit in fits:
        frames = sorted(fit.glob("step_*.dat"))
        assert len(frames) >= 2, (fit, frames)
        for path in frames:
            frame = path.read_bytes()
            assert frame.rstrip().endswith(b"# complete"), path
            outputs[path.relative_to(root)] = frame
    outputs["history"] = [line for line in (root / "history").read_bytes().splitlines()
                          if b"CPU time:" not in line]
    return outputs


def main():
    work, source = (Path(arg) for arg in sys.argv[1:3])
    ranks, command = int(sys.argv[3]), sys.argv[4:]
    with tempfile.TemporaryDirectory(prefix="orientation-cross-call-", dir=work) as tmp:
        experiment = make_experiment(Path(tmp) / "seed", source, command)
        reference_counts = reference_outputs = None
        for mode in ("full", "intra-stats", "stats"):
            case = Path(tmp) / mode
            prepare(case, source, experiment)
            stats = run(case, command, mode)
            assert "MINPACK completed:" in stats, stats
            counts = summarize(stats, ranks, mode)
            outputs = read_outputs(case)
            if reference_outputs is None:
                reference_counts, reference_outputs = counts, outputs
            assert outputs.keys() == reference_outputs.keys(), mode
            for name, contents in outputs.items():
                assert contents == reference_outputs[name], (mode, name)
            assert [len(job) for job in counts] == [len(job) for job in reference_counts], mode
            for job in counts:
                assert all(full + cfonly + spectra == ORIENTATIONS
                           for full, cfonly, spectra, _ in job), (mode, job)
                if mode == "stats":
                    # First residual, Gaussian probe, diffusion probe, then the
                    # ODF-only probe that returns to the base local parameters.
                    assert [row[:3] for row in job[:4]] == [
                        (ORIENTATIONS, 0, 0), (0, ORIENTATIONS, 0),
                        (ORIENTATIONS, 0, 0), (0, 0, ORIENTATIONS),
                    ], job[:4]
                    assert [row[3] for row in job[:4]] == [0, ORIENTATIONS, 0, ORIENTATIONS], job[:4]
                else:
                    assert all(row == (ORIENTATIONS, 0, 0, 0) for row in job), (mode, job)
            print(mode, "initial evaluation counts (full, cfonly, spectra, cross):",
                  [job[:4] for job in counts])


if __name__ == "__main__":
    main()
