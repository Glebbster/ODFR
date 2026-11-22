import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def make_deck(root, source):
    # Start from the smoke deck and keep the legacy positional .par layout.
    # Use a small basis and three angles so this test can run during CTest.
    shutil.copyfile(source / "smoke.esr", root / "smoke.esr")
    odf = (source / "ODFR4F.par").read_text(encoding="ascii").splitlines()
    odf[1] = "0 0.0 80"
    (root / "ODF.par").write_text("\n".join(odf) + "\n", encoding="ascii")
    template = (source / "ODFR4.par").read_text(encoding="ascii").splitlines()
    base = template[:31]
    base[0] = "2 1 0"
    base[1] = "0 2"
    base[2] = "72 2"
    base[4] = "57.0 1.0 0"
    base[9] = "3 0"
    base[12] = "0.001 500"
    base[13] = "6 5 4 3 2"
    spec_lines = []
    for i, odf_type in enumerate((2, 1, 0, 0, 0, 0, 0, 0, 0), 1):
        params = base.copy()
        params[30] = f"0.0 0  {10 * i}.0 0  0.0 0"
        # Use a shared coefficient slot only for the second site of p3.
        params += ["1" if i == 3 else "0", *params[7:31]]
        params[41] = "0.0 0  0.0 0  0.75 10"
        # A flag of 10 keeps a local fixed value in later spectra.
        if i == 2:
            params[27] = "0.05 0  0.06 0  1.0 10"
        if i == 4:
            params[17] = "2.0 10  0.6 0  0.0 0  0.0 0"
        if i == 5:
            params[16] = "0.0 0  25.0 10  1.0 0"
        if i == 6:
            params[2] = "64 2"
        if i == 7:
            params[13] = "8 7 4 3 2"
        if i in (8, 9):
            params[9] = "1 0"
        (root / f"p{i}.par").write_text("\n".join(params) + "\n", encoding="ascii")
        spec_lines.append(f"p{i}.par, r{i}.par, smoke.esr, r{i}.dat, {odf_type}")
    (root / "namODFR.nam").write_text(
        f"1\nODF.par, ODFout.par, {len(spec_lines)}\n"
        + "\n".join(spec_lines) + "\nend\n", encoding="ascii")
    return len(spec_lines)


def run(root, command, mode):
    env = dict(os.environ, ODFR4_ORIENTATION_CACHE=mode)
    result = subprocess.run(command, cwd=root, env=env, timeout=90,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(result.stdout)
    matches = re.findall(r"ORIENTATION_CACHE rank=\d+ full=(\d+) cfonly=(\d+) spectra=(\d+)",
                         result.stdout)
    assert matches, result.stdout
    return tuple(sum(int(row[i]) for row in matches) for i in range(3))


def main():
    work, source = (Path(arg) for arg in sys.argv[1:3])
    command = sys.argv[3:]
    with tempfile.TemporaryDirectory(prefix="orientation-cache-", dir=work) as tmp:
        baseline, cached = Path(tmp) / "full", Path(tmp) / "cached"
        for folder in (baseline, cached):
            folder.mkdir()
            count = make_deck(folder, source)
        full_counts = run(baseline, command, "full")
        cached_counts = run(cached, command, "stats")
        print("FULL, CFONLY, spectra:", full_counts, cached_counts)
        assert full_counts == (50, 0, 0), full_counts
        assert sum(cached_counts) == full_counts[0], cached_counts
        full, cfonly, spectra = cached_counts
        assert full > 0 and cfonly == 0 and spectra > 0, cached_counts
        for name in ("ODFout.par", *(f"r{i}.{ext}" for i in range(1, count + 1)
                                     for ext in ("dat", "par"))):
            assert (baseline / name).read_bytes() == (cached / name).read_bytes(), name


if __name__ == "__main__":
    main()
