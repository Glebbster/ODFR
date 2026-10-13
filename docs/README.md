# Weights and optimization steps / Веса и шаги подбора

[Architecture / Архитектура](ARCHITECTURE.md) · [Manual / Мануал](ODF_manual4_3.pdf)

## Русский

Настройки задаются в `.par` первого спектра, указанного в `namODFR.nam`.
Пример: [s2_2_00.par, строки 1–6](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1-L6).

| Где | Что менять |
|---|---|
| [Строка 5](../examples/ODFR_Tests/sampl2/s2_2_00.par#L5): `4.9  1.00000  4` | Третье число `iwflag`: `0` — без весов; `4` — веса из файла; `5` — веса и деление невязки на уровень шума `rmsn` (второе число этой строки). Остальные значения сохраняют прежние режимы из мануала. |
| [Строка 3](../examples/ODFR_Tests/sampl2/s2_2_00.par#L3): `200  2` | Первое число — количество расчётных точек одного спектра. |
| [Строка 1](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1): `3  0  2` | Первое число `3` выбирает подбор NL2SOL, для которого доступна запись шагов. |

**Веса.** Положите файл `weight` рядом с `namODFR.nam`: одно число на строку, сначала все точки первого спектра, затем второго и т. д. В примере два спектра по 200 точек — нужны 400 чисел. `1` сохраняет вклад точки, `0` исключает её, `2` увеличивает вклад в сумму квадратов ошибки в 4 раза. Файл можно создать в [Weights.py](../tools/README.md). Если файла нет или число значений неверно, программа использует единичные веса и сообщает об этом.

**Шаги подбора.** По умолчанию запись выключена. Для подробного изучения шагов моделирования включите её перед запуском ODFR в том же терминале:

```powershell
$env:ODFR_OPTIMIZATION_STEPS = '1'  # включить; '0' — выключить
```

Начальный и принятые спектры сохраняются в `optimization_steps/fit_*/`. Переключателя в `.par` нет. Linux/macOS: `export ODFR_OPTIMIZATION_STEPS=1` (или `0`).

## English

Edit the `.par` for the first spectrum listed in `namODFR.nam`.
Example: [s2_2_00.par, lines 1–6](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1-L6).

| Location | Setting |
|---|---|
| [Line 5](../examples/ODFR_Tests/sampl2/s2_2_00.par#L5): `4.9  1.00000  4` | Third value, `iwflag`: `0` = no weights; `4` = point weights; `5` = point weights plus division by noise level `rmsn` (second value). Other values retain the modes in the manual. |
| [Line 3](../examples/ODFR_Tests/sampl2/s2_2_00.par#L3): `200  2` | First value = calculated points per spectrum. |
| [Line 1](../examples/ODFR_Tests/sampl2/s2_2_00.par#L1): `3  0  2` | First value `3` selects NL2SOL fitting, for which step recording is available. |

**Weights.** Place `weight` beside `namODFR.nam`, one number per line, all points of the first spectrum followed by the next. This example has two 200-point spectra: 400 values. `1` keeps the contribution, `0` excludes a point, and `2` multiplies its squared-error contribution by 4. Use [Weights.py](../tools/README.md) to create the file. Missing or wrongly sized files produce a diagnostic and unit weights.

**Optimization steps.** Recording is off by default. To examine modeling steps in detail, enable it in the same terminal before running ODFR:

```powershell
$env:ODFR_OPTIMIZATION_STEPS = '1'  # enable; '0' disables
```

Initial and accepted spectra are saved in `optimization_steps/fit_*/`. There is no `.par` switch. Linux/macOS: `export ODFR_OPTIMIZATION_STEPS=1` (or `0`).
