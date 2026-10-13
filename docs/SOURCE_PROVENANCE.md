# Source provenance

- Imported repository: https://github.com/Glebbster/Projects
- Branch snapshot: remote `main`, commit `13f0f3f149e401885ac5fed49c61f4c48cd41bf7`.
- Imported subtree: `Fortran_prog_for_modeling/builds/ODFR4_3_final`.
- Local `main` was `5d8ae70a740b847a8c0725f1fe78f4e458b946ff`; its ODFR subtree matched the remote snapshot. Uncommitted work and other development branches were not imported.
- The imported `VERSION` was `4.3.0`. This standalone distribution is numbered `5.0.0` at the maintainer's request. Scientific computation files match that snapshot except the program name string `ODFR4_3` -> `ODFR5` in `src/ODFR4.f`. Build/output identity and smoke-banner assertions were updated separately.

The functional changelog was checked against the official [SourceForge files](https://sourceforge.net/projects/odfr/files/), downloaded on 2026-10-13:

| Archive | SHA-256 |
|---|---|
| ODFR4_3-Source.zip (serial) | `1174af2e3c763337b144c39688344a27240201097cb12b7d83d0be685d2bca29` |
| ODFRM4_3-Source.zip (MPI) | `6274a702f2567279655ec5a60861c2bfcc0837467ed75489a232f7f1e70904d9` |

The supplied `ODF_manual4_3.pdf` is unchanged, SHA-256 `be137596f1e8d8955df6df00da7a3eba96ea8fc7b80d8597c69d8b8e3bb8334d`. Its old MPI notes describe the earlier implementation.

Examples were copied from the explicitly supplied `ODFR_Tests` directory. Their exact included-file hashes and exclusions are in [`examples/manifest.json`](../examples/manifest.json). The separate imported `tests/` fixtures are the build/regression tests, not substitutes for those examples.

Legacy developer notes under `docs/` retain their historical context. The root README and 5.0.0 changelog are authoritative for release installation and changes.
