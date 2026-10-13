# Changes

[Русский](CHANGELOG_RU.md)

## 5.0.0 — 2026-10-13

First standalone release of the ODFR version maintained in `Glebbster/Projects/main`, imported from `13f0f3f149e401885ac5fed49c61f4c48cd41bf7`. The changes below are cumulative relative to the official **2018 ODFR 4.3 serial/MPI source archives**, not all newly implemented for this release. The internal Projects tag `odfr4-v4.3.0` already contained many of them.

### Functional changes and files

- **Point-weighted fitting:** `src/lfunR.f`. `iwflag=4` multiplies residuals and the explicit Jacobian by coefficients from `weight` / `weight.dat`; `iwflag=5` additionally uses `rmsn`. Requires exactly `ndatot` values; a missing/invalid-size file falls back to unit coefficients with a diagnostic. These are residual multipliers, not direct squared-loss coefficients.
- **Recalculation instead of residual-path cache reuse:** `src/lfunR.f` now sets `icalc=FULL` for each spectrum/site in `iflag=1`. The explicit-Jacobian path (`iflag=2`) still selects `CFONLY` or `FULL` using `modtd`; this is not removal of every cache. Disk storage was already disabled by `DiscTrid=.false.` in the official serial 4.3 code. ODF coefficient changes now invalidate the affected calculations/dependencies.
- **MPI architecture:** `src/ODFR4.f`, `src/momdls.f`, `src/mpi_state.inc`; MPI-aware output in `src/INPNamDR.for`, `src/PARreadN.for`, `src/specread.for`, `src/lfunR.f`, `src/fitNL2.for`, `src/MADR.f`, `src/SpecOut.for`, `src/PARwrN.for`, `src/ODFwrite.for`. Replaces the old `odfrm4_3.for` / `momdlsMPI3#.f` / `CommunucationMPI.f` controller-worker implementation with a common serial/MPI source path. Every rank follows the calculation; orientations are assigned cyclically (`1+rank`, stride `size`) and summed with `MPI_ALLREDUCE`. Normal result writes are restricted to rank 0; single-orientation evaluation runs on rank 0. `compat/serial/mpif.h` and `compat/serial/mpi_stubs.f90` support a build without MPI. No process-count/`NORT` divisibility requirement is encoded in the new loop.
- **Loss diagnostics and best evaluated fit:** `src/lfunR.f`, `src/iterat.inc`, `src/NLSINIT.f`, `src/MADR.f`, `src/fitNL2.for`, `src/SpecOut.for`. Adds `SSE_W`, `SSE_UNW`, reduced losses and per-spectrum contributions; weighted totals are computed after objective normalization. NL2SOL callbacks retain `xbest` at the lowest evaluated weighted SSE, and `fitNL2` restores it after the solver. `SSE_UNW` uses the unweighted residual; `SSE_W` includes enabled weights/normalization/penalties. The new DAT header fields replace the old loss labels; numeric spectrum columns remain. Solver termination and covariance diagnostics still need to be assessed separately.
- **Accepted NL2SOL state recording:** `src/optimization_steps.for`, `src/nl2sol_2015.f90`, `src/fitNL2.for`. Records initial/accepted spectra in fresh `optimization_steps/fit_*` directories on rank 0. Rejected candidates, Jacobian and covariance probes are excluded. Recording is passive and switches off on I/O failure; `ODFR_OPTIMIZATION_STEPS=0` disables it. Tests: `tests/optimization_steps/`.
- **Build and distribution:** `CMakeLists.txt`, `VERSION`, `src/version.inc.in`, `platforms/{windows,linux,macos}/Makefile`, `compat/serial/`, `tests/smoke/`, `tests/hgo_h2_rt/`. Adds a common CMake build, GNU legacy-Fortran flags, serial/MPI targets and regression checks. This release changes `VERSION` to `5.0.0`, the program banner to `ODFR5` (`src/ODFR4.f`), executable output names to `ODFR5_serial` / `ODFR5_mpi`, and the corresponding smoke-test banner assertion. Existing CMake target/option names remain. Includes English/Russian README files, the unchanged manual and the supplied example datasets.

### Fixes and numerical guards

- Initialized ODF state, Wigner coefficients/null flags, loss accumulators, error codes and the potential-term counter; retained orientation storage indices across calls. Removed the callback's use of an uninitialized control flag.
- Preserved the ODF type separately for each task in a multi-task NAM file; checked ODF rank/parameter capacity and handled missing parameter files and malformed spectrum reads more safely.
- Added guards for a one-point field grid, zero orientation/ODF normalization denominators, near-zero noise normalization, empty integration intervals, invalid integration settings/initial NaNs, continued-fraction size/phase/denominator errors, and excessive tridiagonal dimensions. These are specific guards, not a general guarantee against invalid inputs.
- Added modified-Bessel overflow checks. Logarithmic parameter packing now uses `max(parameter,1)` and unpacking uses a shared upper log threshold of 13; these guards can affect extreme/out-of-domain inputs.
- Corrected legacy declarations/COMMON layout and reduced duplicate concurrent writes under MPI.

No robust-covariance prototype or reverted direct-ODF experiment from other branches is included. The current scientific source was imported from main; release preparation itself changes version/build identity, documentation and packaging only. Windows validation: all five CTest cases passed.

- **Weight editor:** `tools/Weights.py` creates regional point weights; [usage](tools/README.md). Its Tkinter interface uses the portable TkAgg backend.
