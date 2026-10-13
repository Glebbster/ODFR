# ODFR 5

[English](README.md) · [Изменения](CHANGELOG_RU.md) · [Релизы](https://github.com/Glebbster/ODFR/releases)

ODFR рассчитывает и подбирает спектры ЭПР спиновых зондов, прежде всего нитроксильных радикалов, в частично ориентированных образцах. Совместный подбор спектров при разных углах образца позволяет исследовать распределение молекулярных ориентаций и вращательную подвижность с помощью моделей MOMD и ODF.

Это развитие [ODFR 4.3](https://sourceforge.net/projects/odfr/), созданной Андреем Х. Воробьёвым. Версия **5.0.0** основана на проверенном снимке исходников `Projects/main`; функциональные изменения и соответствующие файлы перечислены в [CHANGELOG_RU.md](CHANGELOG_RU.md).

[Установка](#установка) · [Мануалы и описание](#мануалы-и-описание) · [Лицензии](#лицензии)

## Установка

### Готовая версия для Windows x64

1. Скачайте **ODFR-5.0.0-windows-x64.zip** из [релиза](https://github.com/Glebbster/ODFR/releases/tag/v5.0.0) и распакуйте.
2. Для `bin/ODFR5_serial.exe` не нужны компилятор, Python или MPI. Библиотеки GNU включены статически.
3. Для `bin/ODFR5_mpi.exe` установите 64-битную среду **Microsoft MPI** (`msmpisetup.exe`) с [сайта Microsoft](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi). SDK для запуска готовой программы не нужен. После установки откройте терминал заново.
4. Создайте отдельную папку расчёта с `namODFR.nam` и указанными в нём входными файлами. Запускайте программу из этой папки:

```powershell
cd C:\work\my-calculation
& 'C:\apps\ODFR-5.0.0\bin\ODFR5_serial.exe'
# Или запуск MPI на двух процессах:
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'C:\apps\ODFR-5.0.0\bin\ODFR5_mpi.exe'
```

Подставьте свои пути распаковки и расчёта. Для одновременных расчётов нужны разные рабочие папки: имена `temp*`, `history` и каталогов записи оптимизации общие. Число процессов выбирайте с учётом свободных ядер и `NORT`; новая реализация MPI не требует кратности числа процессов и `NORT`.

### Сборка из исходников

Нужны **CMake 3.20+, GNU Fortran и инструмент сборки**. Python 3 включает тест наблюдателя оптимизации; самой программе Python не нужен. MPI устанавливается только для параллельной версии. Последовательная версия использует комплектные заглушки MPI; внешние BLAS/LAPACK не требуются.

Скачайте и распакуйте **ODFR-5.0.0-source.zip** либо установите Git и выполните:

```text
git clone https://github.com/Glebbster/ODFR.git
cd ODFR
git checkout v5.0.0
```

Все команды сборки выполняются из папки с `CMakeLists.txt`.

#### Windows: установка пакетов и компиляция

Установите [MSYS2](https://www.msys2.org/), откройте **MSYS2 MINGW64** и обновите пакеты:

```bash
pacman -Syu
```

Если обновление просит закрыть терминал, откройте MINGW64 заново и повторите `pacman -Syu`. Установите 64-битные инструменты:

```bash
pacman -S --needed mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-cmake mingw-w64-x86_64-make mingw-w64-x86_64-python
```

Откройте **PowerShell**. Для стандартного пути установки MSYS2:

```powershell
$env:PATH = "C:\msys64\mingw64\bin;$env:PATH"
cd C:\path\to\ODFR
cmake -S . -B out/windows -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=OFF -DBUILD_TESTING=ON '-DCMAKE_EXE_LINKER_FLAGS=-static' '-DPython3_EXECUTABLE=C:/msys64/mingw64/bin/python.exe'
cmake --build out/windows --parallel 2
ctest --test-dir out/windows --output-on-failure
```

Результат: `out/windows/bin/ODFR5_serial.exe`.

Для MPI установите среду Microsoft MPI по инструкции выше, затем в **MSYS2 MINGW64** добавьте совместимые с MinGW заголовки и библиотеку импорта:

```bash
pacman -S --needed mingw-w64-x86_64-msmpi
```

В **PowerShell** с тем же PATH создайте отдельную сборку:

```powershell
cmake -S . -B out/windows-mpi -G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_SERIAL=ON -DODFR4_BUILD_MPI=ON -DBUILD_TESTING=ON '-DCMAKE_EXE_LINKER_FLAGS=-static' '-DPython3_EXECUTABLE=C:/msys64/mingw64/bin/python.exe' '-DMPIEXEC_EXECUTABLE=C:/Program Files/Microsoft MPI/Bin/mpiexec.exe'
cmake --build out/windows-mpi --parallel 2
ctest --test-dir out/windows-mpi --output-on-failure
```

В `out/windows-mpi/bin/` появятся `ODFR5_serial.exe` и `ODFR5_mpi.exe`. Пакет MSYS2 предоставляет интерфейс для сборки, а отдельно установленная среда Microsoft MPI — `mpiexec.exe` и `msmpi.dll`.

#### Ubuntu / Debian

```bash
sudo apt update
sudo apt install build-essential gfortran cmake python3
cmake -S . -B out/linux -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_MPI=OFF -DBUILD_TESTING=ON
cmake --build out/linux --parallel 2
ctest --test-dir out/linux --output-on-failure
```

Для MPI:

```bash
sudo apt install libopenmpi-dev openmpi-bin
cmake -S . -B out/linux-mpi -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_MPI=ON -DBUILD_TESTING=ON
cmake --build out/linux-mpi --parallel 2
ctest --test-dir out/linux-mpi --output-on-failure
```

Программы находятся в `bin/` выбранной папки сборки. MPI запускается командой `mpiexec -n 2 /path/to/ODFR5_mpi` из папки расчёта.

#### macOS

Установите [Homebrew](https://brew.sh/) по инструкции на его сайте, затем:

```bash
xcode-select --install
brew install gcc cmake make python
```

Проверьте имя компилятора командой `command -v gfortran`. Если установлен только вариант с номером, например `gfortran-15`, подставьте его вместо `gfortran` ниже. Компилятор и MPI должны иметь одинаковую архитектуру.

```bash
cmake -S . -B out/macos -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DODFR4_BUILD_MPI=OFF -DBUILD_TESTING=ON
cmake --build out/macos --parallel 2
ctest --test-dir out/macos --output-on-failure
```

Для MPI сначала выполните `brew install open-mpi`, затем повторите три команды с каталогом `out/macos-mpi` и `-DODFR4_BUILD_MPI=ON`. Программы появятся в `bin/` внутри папки сборки.

**Проверка релиза:** сборка Windows x64 проверена. Команды Linux и macOS описывают конфигурацию сборки, но для этого релиза на этих ОС не выполнялись. Имена параметров `ODFR4_BUILD_*` и целей CMake сохранены для совместимости.

Пакеты: [CMake в MSYS2](https://www.msys2.org/docs/cmake/), [MPI в MSYS2](https://packages.msys2.org/packages/mingw-w64-x86_64-msmpi), [GCC в Homebrew](https://formulae.brew.sh/formula/gcc), [Open MPI в Homebrew](https://formulae.brew.sh/formula/open-mpi).

## Мануалы и описание

- [Архитектура](docs/ARCHITECTURE.md).
- [Мануал](docs/ODF_manual4_3.pdf) — исходный PDF для ODFR 4.3 без изменений. Старые имена программ и описание MPI уточняются этим README и списком изменений.
- [ODFR_Tests](examples/ODFR_Tests) — десять наборов примеров с исходными историческими результатами. Перед запуском скопируйте нужный пример в отдельную папку. Эти результаты не объявляются заново проверенными эталонами ODFR 5.
- [Проверки сборки](tests/README.md) — взвешенный smoke-тест, регрессия HGO, их serial/MPI-варианты, проверка принятых состояний NL2SOL. Перед запуском всего CTest соберите **все цели**; для теста наблюдателя CMake должен найти Python.
- [Веса и шаги подбора](docs/README.md): точные строки PAR и включение записи. [Weights.py](tools/README.md) позволяет задать веса на графике.

Для подробного изучения шагов моделирования можно включить запись перед запуском в том же терминале. По умолчанию она выключена.

```powershell
$env:ODFR_OPTIMIZATION_STEPS = '1'  # включить; '0' — выключить
```

## Лицензии

[GNU GPL v3](LICENSE). Исходные уведомления об авторских правах сохранены. Описание научной модели и сведения для цитирования находятся в неизменённом мануале и на [странице исходного проекта](https://sourceforge.net/projects/odfr/).

[Лицензии библиотек](licenses/README.md).
