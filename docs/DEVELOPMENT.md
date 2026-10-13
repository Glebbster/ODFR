# Разработка

## Сборка

```bash
cmake -S . -B out/local \
  -DODFR4_BUILD_SERIAL=ON \
  -DODFR4_BUILD_MPI=OFF \
  -DBUILD_TESTING=ON
cmake --build out/local --parallel
```

## Тесты

```bash
ctest --test-dir out/local --output-on-failure
```

## Документация

При изменении поведения обновляются:

- `CHANGELOG.md`;
- соответствующий документ в `docs`;
- страница изменённого файла в `docs/reference`.

История разработки связывается только с проверенными Git-коммитами.
