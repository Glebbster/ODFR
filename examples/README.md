# ODFR examples / Примеры ODFR

`ODFR_Tests/` contains the user-supplied ten datasets: `F16nr`, `F322-21br`, `F93-7r`, `k18-6c`, `k18p-1d`, and `sampl1`–`sampl5`. All 456 included files retain their original bytes, including parameter files, spectra, weights and historical outputs. [`manifest.json`](manifest.json) records sizes and SHA-256 hashes.

Excluded: `.DS_Store`, the old `odfr4_3.exe`, and five old macOS `ODFR4_3_weights` executables. Use the new release executables. Original directories were not edited.

Copy a dataset to a separate calculation directory and inspect its `namODFR.nam` before running. Supplied `history` and result files are historical examples, not newly generated ODFR 5 reference outputs. Long example fits are not part of CTest; release checks use the isolated fixtures under `tests/`.

По-русски: сохранены десять переданных наборов и точные байты 456 файлов. Исключены служебный `.DS_Store` и старые программы Windows/macOS; используйте новые бинарники релиза. Перед запуском копируйте пример в отдельную папку и проверьте `namODFR.nam`. Старые результаты сохранены как исторические, а не как заново проверенные эталоны ODFR 5.
