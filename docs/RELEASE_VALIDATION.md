# Release validation — 5.0.0

Verified on 2026-10-13, Windows x64.

- GNU Fortran 15.2.0 (MSYS2 MINGW64); CMake 4.1.1; Python 3.13.0.
- CMake Release configuration, serial + Microsoft MPI, all targets built; `-static` for GNU runtime linkage.
- All five CTest cases passed: `optimization_steps_observer`, `smoke_serial`, `hgo_h2_rt_serial`, `smoke_mpi_2`, `hgo_h2_rt_mpi_2`. Tests ran sequentially, MPI tests used two ranks.
- Separate archive check passed for 14 recorded states: exact reconstruction, bounded-error reconstruction, refusal to overwrite an archive, unchanged source-table hashes.
- PE import inspection: serial depends on `KERNEL32.dll` and `msvcrt.dll`; MPI additionally depends on `msmpi.dll`. No separately installed GNU runtime DLLs are required for these binaries.
- All 456 included example files match their supplied source bytes. The PDF matches the supplied original byte-for-byte.
- Computational source files match the imported main snapshot; the only change within `src/` is the program name string in `ODFR4.f`. The CMake project/output names, VERSION and smoke-test banner assertion were updated.

Build reproduction (from the source root with MinGW tools on PATH and Microsoft MPI installed):

```powershell
cmake -S . -B out/windows-release -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=ON -DBUILD_TESTING=ON '-DCMAKE_EXE_LINKER_FLAGS=-static'
cmake --build out/windows-release --parallel 2
ctest --test-dir out/windows-release --output-on-failure -j 1
```

Python was explicitly selected via `-DPython3_EXECUTABLE=<installed Python 3 path>` in the verified build. The archive check used NumPy from that interpreter. See [test instructions](../tests/README.md) for the separate archive command.

These are regression checks, not performance benchmarks or validation of every model/input combination. The compiler emits warnings from legacy Fortran (including assumed-size array idioms); a warning-free build is not claimed. The ten historical example datasets were packaged but not all rerun. Linux and macOS were not built or run in this release session.

Serial/MPI DAT outputs are byte-identical in both regression cases (`smoke_out.dat` and `final1.dat`), including the HGO rows with no experimental sample (1000 experimental samples and 1024 model points). This comparison concerns the rounded output format, not bitwise identity of internal floating-point calculations.

Microsoft MPI runtime DLL version: `10.1.12498.52`.
