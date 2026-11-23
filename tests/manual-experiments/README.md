# Manual spectrum examples

Each example uses the saved parameters from iteration 6 for one forward
calculation with two sites. Nominal angles are 0, 10, 20, 30, 40, 50, 60,
70, 80, 90, 125, 145, 165 and 180 degrees.

| Directory | Dataset | Point weights | Saved psi offset (degrees) |
| --- | --- | --- | ---: |
| `iso-fl-rt-unweighted-fixed-23p4` | ISO FL at 25 deg C, local C20 = 23.4 at both sites | None | -2.617706 |
| `lc-fl-weighted-psi-it6` | LC FL at 25 deg C, weighted psi fit | `weight` | -5.329084 |
| `lc-p-weighted-psi-it6` | LC P at 25 deg C, weighted psi fit | `weight` | -0.374605 |

Copy one complete example directory to a fresh writable location. Keep
`namODFR.nam`, `data/`, `params/`, `odf.par`, `outputs/` and, where present,
`weight` together. The source data are the denoised experimental ESR files
used for the fit.

On Windows, open a terminal in the copy and use the absolute path to your
executable. For example:

```powershell
Set-Location 'C:\work\iso-fl-rt-unweighted-fixed-23p4'
& 'C:\apps\ODFR\bin\ODFR5_serial.exe'
```

For an MPI build, use its matching MPI launcher. With Microsoft MPI:

```powershell
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'C:\apps\ODFR\bin\ODFR5_mpi.exe'
```

On macOS or Linux, use the executable path without the `.exe` extension.
Run from the same example directory.

The program reads `namODFR.nam` from the current directory. The
`Kcontrol=2` in the first PAR selects a single forward calculation with
zero fitted variables. The input fit flags are retained. The program reads
the experimental spectra for the loss report. The LC examples use the
included point weights and the experimental noise setting in the PAR
files. The `weight` file contains 1,024 entries per spectrum in NAM order.
The ISO example has point weighting and the ODF penalty disabled.

A successful run marks the simulation `Finished` in `history`, reports
`N_TOTAL=14336` and `NPRM=0`, and produces 14 DAT and 14 PAR files in
`outputs/`. It also writes `outputs/odf.par`,
`G1outputs/odf.par` and `G2outputs/odf.par`. Use the loss report in
`history` or the DAT headers to compare results.

All three examples use a model field grid of 3302 to 3462 G with 1,024
points per spectrum. The values below are approximate.

| Dataset | Expected total SSE_W |
| --- | ---: |
| ISO FL | 0.0000944297 |
| LC FL | 0.0461690 |
| LC P | 0.0300660 |

For ISO FL, SSE_W equals SSE_UNW. For each LC example, check that the
log also reports `Weights OK, ndatot = 14336`.
