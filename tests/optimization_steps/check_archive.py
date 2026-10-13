"""Check exact archives, bounded integer deltas, and unchanged source tables."""
import hashlib
import importlib.util
from pathlib import Path
import sys
import tempfile

import numpy as np

script = Path(__file__).resolve().parents[2] / "tools" / "optimization_animation.py"
spec = importlib.util.spec_from_file_location("animation", script)
animation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(animation)
source = Path(sys.argv[1]).resolve()
root = Path(tempfile.mkdtemp(prefix="archive_check_", dir=Path(sys.argv[2]).resolve()))
files = sorted(source.glob("*.dat"))
before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}

for precision, name in ((None, "exact.zip"), (1e-6, "compact.zip")):
    output = root / name
    animation.pack(source, output, precision)
    original = animation.History(source)
    restored = animation.History(output)
    assert original.count == restored.count
    assert np.array_equal(original.experiment, restored.experiment)
    peaks = {int(i): v for i, v in restored.manifest["spectrum_peak_amplitude"].items()}
    measured = 0.0
    for (left_meta, left), (right_meta, right) in zip(original, restored):
        assert left_meta == right_meta
        assert np.array_equal(left[:, :2], right[:, :2])
        if precision is None:
            assert np.array_equal(left, right)
        else:
            bounds = np.array([(peaks[int(i)] or 1)*precision/2 for i in left[:, 0]])
            error = np.abs(left[:, 2] - right[:, 2])
            assert np.all(error <= bounds * (1+1e-8)), float(np.max(error/bounds))
            measured = max(measured, float(np.max(error)))
    print(f"PASS: {name}, {original.count} frames, max error {measured:g}")
    try:
        animation.pack(source, output, precision)
    except FileExistsError:
        pass
    else:
        raise AssertionError("Existing archive was overwritten")

assert before == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
print(f"PASS: source table hashes unchanged; artifacts: {root}")
