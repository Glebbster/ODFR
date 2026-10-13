# Архитектура

Основная последовательность вызовов:

```text
ODFR4
├─ INPNamDR
├─ NLSINIT
├─ PARreadN / ODFread
├─ specread / choiceExp
├─ прямой расчёт: lfun
├─ подбор: fitl или fitNL2
│  └─ lfun
│     └─ momdls
│        └─ eprls
└─ SpecOut / PARWrite / ODFwrite
```

Последовательная и MPI-сборки используют один набор файлов из `src`. Последовательная цель подключает заглушки из `compat/serial`.
