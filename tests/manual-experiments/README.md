# Membrane EPR examples

These examples simulate how electron paramagnetic resonance (EPR) spectra of
membrane samples vary with orientation. The datasets come from an unpublished
manuscript. The DOI will be added after publication.

| Example directory | Membrane preparation | Spectral point weights |
| --- | --- | --- |
| `lc-fl-weighted` | Filtration from the liquid-crystalline phase (LC-FL) | Included |
| `lc-p-weighted` | Liquid-crystalline membrane prepared under pressure (LC-P) | Included |
| `iso-fl-rt-unweighted-fixed-23p4` | Filtration from the isotropic phase (ISO-FL) | None |

Each run calculates spectra for 14 sample orientations without fitting.
The macroscopic orientational distribution in `odf.par` uses terms through
order 8.

## Run an example

1. Copy the entire example directory to a writable location.
2. Open PowerShell and run the ODFR executable from inside that copy, where `namODFR.nam`
   is located. Adjust both paths for your installation and example.

   ```powershell
   Set-Location 'C:\work\lc-fl-weighted'
   & 'C:\apps\ODFR\bin\ODFR5_serial.exe'
   ```

For the Windows MPI build, use its matching launcher from the example folder.
With Microsoft MPI:

```powershell
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'C:\apps\ODFR\bin\ODFR5_mpi.exe'
```

On macOS or Linux, run the built executable from the same directory using
its absolute path without the `.exe` extension.

## View the spectra

Check `history` for `Finished`. Open a text file such as `outputs/000.dat`
and plot the two spectra as separate XY series using a spreadsheet or a
plotting tool. The DAT filenames indicate the nominal sample angles in
degrees.

Use `Hexp` and `Sexp` for the experimental curve and `Hmod` and `Smod` for
the total model curve. Each curve has its own magnetic-field column.
