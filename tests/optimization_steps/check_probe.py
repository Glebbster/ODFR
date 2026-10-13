"""Run the analytic NL2SOL probe with recording on/off and with an I/O failure.

Usage: python check_probe.py path/to/optimization_steps_probe.exe output_parent
Only creates a new temporary directory within output_parent; never edits inputs.
"""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def read_step(path):
    lines = path.read_text().splitlines()
    assert lines[-1] == "# complete", path
    headers = dict(line[2:].split(maxsplit=1) for line in lines[:-1]
                   if line.startswith("# "))
    rows = [[float(v) for v in line.split()] for line in lines
            if not line.startswith("#")]
    return headers, rows


def main():
    exe = str(Path(sys.argv[1]).resolve())
    parent = Path(sys.argv[2]).resolve()
    parent.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="probe_", dir=parent))
    for mode in ("off", "on", "failure"):
        run = root / mode
        run.mkdir()
        if mode == "failure":
            (run / "optimization_steps").write_text("block directory creation")
        env = dict(os.environ, ODFR_OPTIMIZATION_STEPS="0" if mode == "off" else "1")
        result = subprocess.run([exe], cwd=run, env=env, capture_output=True, text=True, check=True)
        (run / "stdout.log").write_text(result.stdout)
        if mode == "on":
            assert "disabled" not in result.stdout, result.stdout
        if mode == "failure":
            assert "disabled" in result.stdout, result.stdout
    for name in ("evaluations.dat", "solution.dat"):
        expected = (root / "off" / name).read_bytes()
        for mode in ("on", "failure"):
            assert (root / mode / name).read_bytes() == expected, (mode, name)
    assert not (root / "off" / "optimization_steps").exists()
    recorded = root / "on" / "optimization_steps" / "fit_000001"
    frames = sorted(recorded.glob("step_*.dat"))
    assert len(frames) > 5
    calls = [[float(v) for v in row.split()] for row in
             (root / "on" / "evaluations.dat").read_text().splitlines()]
    ordinary = {}
    for row in calls:
        if row[1] == 12 or (row[1] == 1 and row[3] == 0):
            # Initial finite differences retain IV(1)=12 and NF=1.
            ordinary.setdefault(int(row[0]), row[4:6])
    previous = math.inf
    for path in frames:
        header, rows = read_step(path)
        x, y = map(float, header["x_optimizer"].split())
        assert [x, y] == ordinary[int(header["evaluation"])]
        model = [-10 * (y - x*x), x - 1, x + y]
        assert [r[0] for r in rows] == [1, 1, 2]
        assert [r[1] for r in rows] == [1, 2, 3]
        for row, value in zip(rows, model):
            assert math.isclose(row[2], value, rel_tol=1e-13, abs_tol=1e-14)
        loss = float(header["SSE_W"])
        assert math.isclose(loss, sum(v*v for v in model[:2]), rel_tol=1e-13, abs_tol=1e-24)
        assert math.isclose(float(header["SSE_UNW"]), sum(v*v for v in model), rel_tol=1e-13)
        assert loss <= previous, (path, loss, previous)
        previous = loss
    assert previous < 1e-14
    assert len(frames) < len(ordinary), "fixture must exercise rejected candidates"
    solution = (root / "on" / "solution.dat").read_text().splitlines()
    assert int(solution[0].split()[4]) > 0, "fixture must exercise covariance probes"
    assert list(map(float, read_step(frames[-1])[0]["x_optimizer"].split())) == list(map(float, solution[1].split()))
    # A second fit in the same directory must leave the first history intact.
    before = {p.name: p.read_bytes() for p in recorded.iterdir()}
    subprocess.run([exe], cwd=root / "on", env=dict(os.environ, ODFR_OPTIMIZATION_STEPS="1"),
                   check=True, capture_output=True)
    assert (recorded.parent / "fit_000002").is_dir()
    assert before == {p.name: p.read_bytes() for p in recorded.iterdir()}
    print(f"PASS: {len(frames)} states, {len(ordinary)} candidate evaluations, {len(calls)} total calls")
    print(f"PASS: identical evaluation trajectory/results with recorder on, off, and I/O failure; {root}")


if __name__ == "__main__":
    main()
