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

## Examples / Примеры

The ten supplied datasets are in [`examples/ODFR_Tests`](../examples/ODFR_Tests). Their historical outputs are not newly certified reference results.

По-русски: соберите все цели и запустите CTest командой выше. Проверяются smoke и HGO для serial/MPI и запись принятых состояний NL2SOL. Нужен Python 3; расчёты выполняются на копиях входов внутри папки сборки.
