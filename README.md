# ODFR 5

[Русский](README_RU.md) · [Changes](CHANGELOG.md) · [Releases](https://github.com/Glebbster/ODFR/releases)

ODFR simulates and fits slow-motion EPR spectra of spin probes, especially nitroxide radicals, in partially oriented samples. It jointly fits spectra recorded at different sample angles to study molecular orientation distributions and rotational mobility, using MOMD and ODF models.

This is a maintained continuation of [ODFR 4.3](https://sourceforge.net/projects/odfr/), originally developed by Andrey Kh. Vorobiev. Version **5.0.0** packages the verified `Projects/main` source snapshot; [CHANGELOG.md](CHANGELOG.md) identifies the functional changes and their source files.

## Download and run on Windows x64

1. Download and extract **ODFR-5.0.0-windows-x64.zip** from [Releases](https://github.com/Glebbster/ODFR/releases/tag/v5.0.0).
2. `bin/ODFR5_serial.exe` needs no compiler, Python or MPI installation. The GNU runtime is linked statically.
3. For `bin/ODFR5_mpi.exe`, install the x64 **Microsoft MPI runtime** (`msmpisetup.exe`) from [Microsoft](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi). The SDK is not needed to run this binary. Reopen the terminal after installation.
4. Prepare a separate calculation directory containing `namODFR.nam` and the input files it names. Run from that directory, using the executable's full path:

```powershell
cd C:\work\my-calculation
& 'C:\apps\ODFR-5.0.0\bin\ODFR5_serial.exe'
# Or use 2 MPI processes:
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'C:\apps\ODFR-5.0.0\bin\ODFR5_mpi.exe'
```

Adjust both paths to your extraction/calculation directories. Each simultaneous calculation needs its own working directory: `temp*`, `history` and optimization recordings use shared names. Choose a process count suited to the available cores and `NORT`; this MPI implementation does not require divisibility between the process count and `NORT`.

## Build from source

Required: **CMake 3.20+, GNU Fortran and a build tool**. Python 3 enables the optimizer observer test; the program itself does not need Python. MPI is optional. Serial builds use included MPI stubs and need no MPI installation or external BLAS/LAPACK.

Download **ODFR-5.0.0-source.zip** and extract it, or install Git and run:

```text
git clone https://github.com/Glebbster/ODFR.git
cd ODFR
git checkout v5.0.0
```

Run the following commands from the directory containing `CMakeLists.txt`.

### Windows: install packages, then compile

Install [MSYS2](https://www.msys2.org/), open **MSYS2 MINGW64**, and update it:

```bash
pacman -Syu
```

If prompted to close the terminal, reopen MINGW64 and repeat `pacman -Syu`. Install the native 64-bit toolchain:

```bash
pacman -S --needed mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-cmake mingw-w64-x86_64-make mingw-w64-x86_64-python
```

Open PowerShell. With the default MSYS2 installation location:

```powershell
$env:PATH = "C:\msys64\mingw64\bin;$env:PATH"
cd C:\path\to\ODFR
cmake -S . -B out/windows -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=OFF -DBUILD_TESTING=ON '-DCMAKE_EXE_LINKER_FLAGS=-static' '-DPython3_EXECUTABLE=C:/msys64/mingw64/bin/python.exe'
cmake --build out/windows --parallel 2
ctest --test-dir out/windows --output-on-failure
```

Result: `out/windows/bin/ODFR5_serial.exe`.

For MPI, install the Microsoft MPI runtime as above and add the MinGW-compatible MPI headers/import library in **MSYS2 MINGW64**:

```bash
pacman -S --needed mingw-w64-x86_64-msmpi
```

Then configure a separate build in **PowerShell** with the same PATH:

```powershell
cmake -S . -B out/windows-mpi -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=ON -DBUILD_TESTING=ON '-DCMAKE_EXE_LINKER_FLAGS=-static' '-DPython3_EXECUTABLE=C:/msys64/mingw64/bin/python.exe' '-DMPIEXEC_EXECUTABLE=C:/Program Files/Microsoft MPI/Bin/mpiexec.exe'
cmake --build out/windows-mpi --parallel 2
ctest --test-dir out/windows-mpi --output-on-failure
```

Result: `out/windows-mpi/bin/ODFR5_serial.exe` and `ODFR5_mpi.exe`. The MSYS2 MPI package supplies the build interface; the separately installed Microsoft runtime supplies `mpiexec.exe` and `msmpi.dll`.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install build-essential gfortran cmake python3
cmake -S . -B out/linux -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_MPI=OFF -DBUILD_TESTING=ON
cmake --build out/linux --parallel 2
ctest --test-dir out/linux --output-on-failure
```

For MPI:

```bash
sudo apt install libopenmpi-dev openmpi-bin
cmake -S . -B out/linux-mpi -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_MPI=ON -DBUILD_TESTING=ON
cmake --build out/linux-mpi --parallel 2
ctest --test-dir out/linux-mpi --output-on-failure
```

Executables are in the selected build directory's `bin/`. Run MPI as `mpiexec -n 2 /path/to/ODFR5_mpi` from the calculation directory.

### macOS

Install [Homebrew](https://brew.sh/) using its official instructions, then:

```bash
xcode-select --install
brew install gcc cmake make python
```

Use `command -v gfortran` to check the compiler name; some installations expose a versioned name such as `gfortran-15`, which must then replace `gfortran` below. Keep the compiler and MPI on the same architecture.

```bash
cmake -S . -B out/macos -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_MPI=OFF -DBUILD_TESTING=ON
cmake --build out/macos --parallel 2
ctest --test-dir out/macos --output-on-failure
```

For MPI, first run `brew install open-mpi`, then repeat the three commands using `out/macos-mpi` and `-DODFR4_BUILD_MPI=ON`. Executables are in `bin/` inside the build directory.

**Release verification:** Windows x64 is tested for this release. Linux and macOS commands describe the supported build configuration; they have not been executed for this release. The existing `ODFR4_BUILD_*` option and CMake target names are retained for compatibility.

## Tests, manual and weights

- [Original manual](docs/ODF_manual4_3.pdf): the supplied ODFR 4.3 PDF, unchanged. Its old executable names and MPI implementation notes must be read alongside this README and the changelog.
- [ODFR_Tests](examples/ODFR_Tests): ten example datasets, including their supplied historical outputs. Copy a case to a separate directory before running it. These are examples, not newly certified reference results.
- [Build tests](tests/README.md): weighted smoke and HGO regression tests, serial and two-process MPI variants, optimizer observer checks. Build **all targets** before running all CTest tests; Python must be found for the observer test.
- Point weights: `iwflag=4` applies point coefficients, `iwflag=5` also applies the existing `rmsn` normalization. `weight` (fallback: `weight.dat`) contains exactly `ndatot` coefficients, one per line, in the concatenated fitted-spectrum grid order. Coefficients multiply residuals, so their contribution to squared loss is squared; zero excludes a point. Missing or wrongly sized files fall back to unit coefficients and print a diagnostic. Existing ODF penalties remain part of the objective.
- [Weights and optimization steps](docs/README.md): exact PAR lines and recording settings. [Weights.py](tools/README.md) creates a weight file interactively.

## License and origin

[GNU GPL v3](LICENSE). Original copyright notices are preserved. The scientific model and citation information are described in the unchanged manual and [original project](https://sourceforge.net/projects/odfr/).

Package references: [MSYS2 CMake](https://www.msys2.org/docs/cmake/), [MSYS2 MPI package](https://packages.msys2.org/packages/mingw-w64-x86_64-msmpi), [Homebrew GCC](https://formulae.brew.sh/formula/gcc), [Homebrew Open MPI](https://formulae.brew.sh/formula/open-mpi).
