# Шрифты IBM Plex

Шрифты используются графиками NeoQC и встраиваются в HTML-отчёт, чтобы текст
выглядел одинаково в любой системе и отчёт работал без интернета.

| Файлы | Потребитель |
|---|---|
| `IBMPlexSans-Regular.ttf`, `IBMPlexSans-Bold.ttf`, `IBMPlexMono-Regular.ttf` | Matplotlib (`scripts/plot_style.py`) — расчёт размеров текста на графиках |
| `IBMPlex{Sans,Mono}-<вес>.{latin,cyrillic}.woff2` | HTML-отчёт (`scripts/qc_report.py`) — встраиваются как `@font-face` |
| `unicode_ranges.json` | диапазоны Unicode подмножеств `latin` и `cyrillic` для `@font-face` |
| `OFL.txt` | лицензия |

Веса: IBM Plex Sans — 400, 500, 600, 700; IBM Plex Mono — 400, 500, 600.

Источник — Google Fonts (IBM Plex, © 2017 IBM Corp.), файлы не изменялись.
Текст лицензии взят из официального пакета `@ibm/plex-sans`.

Лицензия — SIL Open Font License 1.1 (`OFL.txt`): шрифты можно свободно
использовать, встраивать и распространять вместе с NeoQC; продавать сами
шрифты отдельно нельзя; «Plex» — зарезервированное имя шрифта, изменённые
версии должны называться иначе.
