# Build and regression tests / Проверки сборки

Build all targets, then run `ctest --test-dir out/windows-mpi --output-on-failure` (substitute your build directory). With Python 3 and both serial and MPI enabled, the default suite has nine tests. The optional stability suite adds two tests.


| Test | Checks |
|---|---|
| `smoke_serial` | 512 supplied point weights, version banner, required outputs and loss-report fields |
| `hgo_h2_rt_serial` | 1024-point HGO case, output contents/size, weighted SSE in the existing `1.xE-02` range |
| `smoke_mpi_2` | Same smoke assertions using two MPI processes |
| `hgo_h2_rt_mpi_2` | Same HGO assertions using two MPI processes |
| `optimization_steps_observer` | Analytic NL2SOL case: accepted states, rollback, rejection/covariance exclusion, recording off by default (unset/0/empty/invalid), enabled by 1, identical trajectory/results in all modes and on I/O failure, preservation of an existing session |
| `user_batch_serial`, `user_batch_mpi_2` | Match a 14-spectrum batch and a second job in the same process to isolated references |
| `user_output_paths_serial`, `user_output_paths_mpi_2` | Create ODF/G1/G2 folders, repeat output, and diagnose blocked creation and opening |
| `user_series_serial`, `user_series_mpi_2` (optional) | Run ten jobs in one process with alternating spectrum counts and CG tolerances; compare all 15 spectra with isolated runs |

CTest copies input datasets into the build directory. Do not point test runners at a live calculation. `ctest -L serial` and `ctest -L mpi` omit the observer test; running all tests avoids this omission. MPI tests need `mpiexec` and a working runtime. These tests verify specific behavior, not the full physical model or performance.

## Windows and repeated runs / Windows и повторы

From the release repository root, with CMake/CTest 3.26+, GNU Fortran, Python 3,
MinGW Make and Microsoft MPI installed:

```powershell
.\platforms\windows\verify.ps1 -BuildDir out/windows-check -Repeat 2 -Stability
```

Use `-ToolchainBin` to select your MSYS2 compiler directory, and
`-PythonExecutable` if Python lives outside that directory. The script builds
both `ODFR5_serial.exe` and `ODFR5_mpi.exe`. Every test must pass all requested
repeats. The first failure fails the run. `-Stability` enables
`ODFR4_TEST_STABILITY`; after configuring and building, use `ctest -L stability`
with the same build directory to run just the two ten-job tests.

GitHub Actions runs the default suite for PRs on Windows x64. Once the workflow
is on the default branch, its manual run enables the stability suite and lets
you select the repeat count. Longer runs for memory trends and explicit
interruption/recovery tests are separate work.

На Windows скрипт собирает обе версии ODFR 5 и запускает тесты.
`-Repeat` задаёт число последовательных проходов каждого теста.
`-Stability` добавляет два сценария по десять заданий в одном процессе.
Обычный PR запускает девять тестов. Ручной запуск в Actions после принятия
workflow добавляет сценарии стабильности. Для длительной проверки памяти
и восстановления после прерывания нужны отдельные сценарии.

## Examples / Примеры

The ten supplied datasets are in [`examples/ODFR_Tests`](../examples/ODFR_Tests). Their historical outputs are not newly certified reference results.

По-русски: соберите все цели и запустите CTest командой выше. Нужен Python 3. Расчёты выполняются на копиях входов внутри папки сборки.
