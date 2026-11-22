[CmdletBinding()]
param(
    [Parameter(Mandatory)][string]$BuildDir,
    [string]$ToolchainBin = 'C:\msys64\mingw64\bin',
    [string]$PythonExecutable,
    [switch]$Stability,
    [ValidateRange(1, 50)][int]$Repeat = 1
)

$ErrorActionPreference = 'Stop'
$source = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..\..')).Path
$bin = (Resolve-Path -LiteralPath $ToolchainBin).Path
$mpiexec = Join-Path $env:ProgramFiles 'Microsoft MPI\Bin\mpiexec.exe'
$compiler = Join-Path $bin 'gfortran.exe'
$mpiCompiler = Join-Path $bin 'mpif90.exe'
if (-not $PythonExecutable) { $PythonExecutable = Join-Path $bin 'python.exe' }

foreach ($tool in @($compiler, $mpiCompiler, $mpiexec, $PythonExecutable)) {
    if (-not (Test-Path -LiteralPath $tool -PathType Leaf)) {
        throw "Required Windows tool not found: $tool"
    }
}

$stabilityTests = if ($Stability) { 'ON' } else { 'OFF' }
$previousPath = $env:Path
try {
    $env:Path = "$bin;$previousPath"
    $configure = @(
        '-S', $source, '-B', $BuildDir, '-G', 'MinGW Makefiles'
        '-DCMAKE_BUILD_TYPE=Release', "-DCMAKE_Fortran_COMPILER=$compiler"
        '-DCMAKE_EXE_LINKER_FLAGS=-static', '-DODFR4_BUILD_SERIAL=ON'
        '-DODFR4_BUILD_MPI=ON', '-DBUILD_TESTING=ON'
        "-DMPI_Fortran_COMPILER=$mpiCompiler", "-DMPIEXEC_EXECUTABLE=$mpiexec"
        "-DPython3_EXECUTABLE=$PythonExecutable", "-DODFR4_TEST_STABILITY=$stabilityTests"
    )
    cmake @configure
    if ($LASTEXITCODE) { throw "CMake configure failed ($LASTEXITCODE)" }

    cmake --build $BuildDir --parallel 4
    if ($LASTEXITCODE) { throw "Windows build failed ($LASTEXITCODE)" }

    ctest --test-dir $BuildDir --parallel 2 --output-on-failure --no-tests=error --repeat "until-fail:${Repeat}"
    if ($LASTEXITCODE) { throw "Windows CTest failed ($LASTEXITCODE)" }
} finally {
    $env:Path = $previousPath
}
