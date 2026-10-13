# Weights.py

## Русский

Оконная программа для задания весов отдельным участкам спектра.
Нужен Python 3 с Tkinter (в установщике Python для Windows включите Tcl/Tk).

```text
python -m pip install numpy pandas matplotlib
```

Перейдите в папку расчёта и запустите скрипт, указав свой путь:

```text
cd C:\work\my-calculation
python C:\path\to\ODFR\tools\Weights.py
```

1. «Считать спектр»: откройте выходной `.dat` ODFR с колонками `Hexp`, `Sexp`, `Hmod`, `Smod`, затем «Отрисовать».
2. Добавьте границы; перетаскивайте их левой кнопкой мыши, удаляйте правой.
3. Введите коэффициенты областей слева направо через запятую. Для двух границ нужны три коэффициента, например `1,0,2`.
4. «N точек» — первое число [строки 3 `.par`](../examples/ODFR_Tests/sampl2/s2_2_00.par#L3). «Копий» повторяет эти веса для нескольких спектров с одинаковой сеткой и областями.
5. «Сгенерировать weight» сохраняет **и заменяет** `weight` в папке запуска. Коэффициенты округляются до одного знака после запятой.

Скрипт строит равномерную сетку от минимального до максимального `Hexp`: этот диапазон должен соответствовать расчётной сетке. Для примера sampl2: N=200, копий=2. [Включение весов в ODFR](../docs/README.md).

Linux: команда запуска `python3`, при отсутствии Tkinter установите `python3-tk`. На macOS используйте Python с Tkinter; проверьте его командой `python3 -m tkinter`.

## English

A windowed tool for assigning weights to spectral regions. Requires Python 3 with Tkinter and `numpy pandas matplotlib` (installation command above). Run the script from your calculation directory.

Open an ODFR output `.dat`, plot it, add/drag region borders, and enter one coefficient per region from left to right. Set **N** to the first value of [PAR line 3](../examples/ODFR_Tests/sampl2/s2_2_00.par#L3); **copies** repeats the pattern for spectra with the same grid/regions. Generate `weight`: it **overwrites the file in the current directory**, with one decimal place. The grid spans min/max `Hexp`, which must match your calculation range. See [ODFR settings](../docs/README.md).
