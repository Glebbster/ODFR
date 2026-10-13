# Build and regression tests / Проверки сборки

Build all targets, then run `ctest --test-dir out/windows-mpi --output-on-failure` (substitute your build directory). CMake must find Python 3 to enable the observer test. The expected full serial+MPI suite has five tests:

| Test | Checks |
|---|---|
| `smoke_serial` | 512 supplied point weights, version banner, required outputs and loss-report fields |
| `hgo_h2_rt_serial` | 1024-point HGO case, output contents/size, weighted SSE in the existing `1.xE-02` range |
| `smoke_mpi_2` | Same smoke assertions using two MPI processes |
| `hgo_h2_rt_mpi_2` | Same HGO assertions using two MPI processes |
| `optimization_steps_observer` | Analytic NL2SOL case: accepted states, rollback, rejection/covariance exclusion, identical trajectory/results with recording on/off/I/O failure, preservation of an existing session |

CTest copies input datasets into the build directory. Do not point test runners at a live calculation. `ctest -L serial` and `ctest -L mpi` omit the observer test; running all tests avoids this omission. MPI tests need `mpiexec` and a working runtime. These tests verify specific behavior, not the full physical model or performance.

## Archive test / Проверка архива

`tests/optimization_steps/check_archive.py` is included but is **not registered in CTest**. It checks lossless reconstruction, bounded-error integer deltas, overwrite rejection and unchanged input hashes. Install NumPy; Matplotlib and Pillow are needed for the separate rendering commands:

```text
python -m venv .venv
```

Activate with `.venv\Scripts\Activate.ps1` in PowerShell, or `source .venv/bin/activate` on Linux/macOS, then:

```text
python -m pip install numpy matplotlib pillow
```

On Windows, after building all targets (use a fresh fixture directory):

```powershell
cmake -E make_directory out/archive-fixture
cmake -E chdir out/archive-fixture ../windows-mpi/bin/optimization_steps_probe.exe
python tests/optimization_steps/check_archive.py out/archive-fixture/optimization_steps/fit_000001 out
```

On Linux/macOS replace `../windows-mpi/bin/optimization_steps_probe.exe` with the corresponding build path and omit `.exe`. Ensure `ODFR_OPTIMIZATION_STEPS` is unset or equals `1` when generating the fixture.

## Historical examples / Исторические примеры

The user-supplied `ODFR_Tests` datasets are in [`examples/ODFR_Tests`](../examples/ODFR_Tests); [`examples/README.md`](../examples/README.md) describes provenance and exclusions. They are not registered as automated CTest cases and were not all rerun for release validation.

По-русски: CTest включает smoke и HGO для serial/MPI, а при найденном Python — проверку наблюдателя NL2SOL (всего пять тестов). Проверка архивирования запускается отдельно командами выше и требует NumPy. Исходные наборы `ODFR_Tests` сохранены как примеры с историческими результатами; они не объявляются пройденными автоматическими тестами нового релиза.
