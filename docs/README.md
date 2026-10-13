# Weights and optimization steps / Веса и шаги подбора

[Architecture / Архитектура](ARCHITECTURE.md) · [Manual / Мануал](ODF_manual4_3.pdf)

## Русский

Настройки задаются в `.par` первого спектра, указанного в `namODFR.nam`.
Пример: [s2_2_00.par, строки 1–6](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1-L6).

| Где | Что менять |
|---|---|
| [Строка 5](../examples/ODFR_Tests/sampl2/s2_2_00.par#L5): `4.9  1.00000  4` | Третье число `iwflag`: `0` — без весов; `4` — веса из файла; `5` — веса и деление невязки на уровень шума `rmsn` (второе число этой строки). Остальные значения сохраняют прежние режимы из мануала. |
| [Строка 3](../examples/ODFR_Tests/sampl2/s2_2_00.par#L3): `200  2` | Первое число — количество расчётных точек одного спектра. |
| [Строка 1](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1): `3  0  2` | Первое число `3` выбирает подбор NL2SOL, при котором записываются шаги. |

**Веса.** Положите файл `weight` рядом с `namODFR.nam`: одно число на строку, сначала все точки первого спектра, затем второго и т. д. В примере два спектра по 200 точек — нужны 400 чисел. `1` сохраняет вклад точки, `0` исключает её, `2` увеличивает вклад в сумму квадратов ошибки в 4 раза. Файл можно создать в [Weights.py](../tools/README.md). Если файла нет или число значений неверно, программа использует единичные веса и сообщает об этом.

**Шаги подбора.** Запись включена по умолчанию для NL2SOL. В `optimization_steps/fit_000001/` сохраняются эксперимент (`experiment.dat`), начальный и принятые расчётные спектры (`step_*.dat`); неудачные пробные шаги не записываются. Следующий запуск создаёт новый номер папки. Отдельного переключателя в `.par` нет. Перед запуском ODFR в том же терминале:

```powershell
$env:ODFR_OPTIMIZATION_STEPS = '1'  # включить; '0' — выключить
```

Linux/macOS: `export ODFR_OPTIMIZATION_STEPS=1` (или `0`).

## English

Edit the `.par` for the first spectrum listed in `namODFR.nam`.
Example: [s2_2_00.par, lines 1–6](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1-L6).

| Location | Setting |
|---|---|
| [Line 5](../examples/ODFR_Tests/sampl2/s2_2_00.par#L5): `4.9  1.00000  4` | Third value, `iwflag`: `0` = no weights; `4` = point weights; `5` = point weights plus division by noise level `rmsn` (second value). Other values retain the modes in the manual. |
| [Line 3](../examples/ODFR_Tests/sampl2/s2_2_00.par#L3): `200  2` | First value = calculated points per spectrum. |
| [Line 1](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1): `3  0  2` | First value `3` selects NL2SOL fitting, which records optimization steps. |

**Weights.** Place `weight` beside `namODFR.nam`, one number per line, all points of the first spectrum followed by the next. This example has two 200-point spectra: 400 values. `1` keeps the contribution, `0` excludes a point, and `2` multiplies its squared-error contribution by 4. Use [Weights.py](../tools/README.md) to create the file. Missing or wrongly sized files produce a diagnostic and unit weights.

**Optimization steps.** Recording is on by default for NL2SOL. `optimization_steps/fit_000001/` contains `experiment.dat` and initial/accepted calculated spectra in `step_*.dat`; rejected trials are omitted. Later runs use a new folder number. There is no `.par` switch. Before running ODFR in the same terminal, use `$env:ODFR_OPTIMIZATION_STEPS = '1'` in PowerShell, or `export ODFR_OPTIMIZATION_STEPS=1` on Linux/macOS; use `0` to disable recording.
