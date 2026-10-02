# Контракт графиков NeoQC

Графики строятся одним процессом Python. C++-ядро записывает стабильные TSV,
`scripts/plot_results.py` проверяет их и строит все доступные метрики, затем
рассчитывает `qc_evaluation.json` и собирает HTML-отчёт.

## Результаты

При запуске с `--plot` в `<result>/plots/` создаются:

| Файл | Назначение |
|---|---|
| `<график>_<рид>.svg` | основной график на английском, векторный; текст остаётся текстом |
| `<график>_<рид>.ru.svg` | тот же график с русскими подписями (для HTML-отчёта) |
| `<график>_<рид>.png` | растровая копия на английском: 2400 × 1350 px, 300 dpi, для экспорта и печати |
| `plots_manifest.json` | машинно-читаемый контракт для отчёта и внешних потребителей |

Рядом с каталогом графиков создаются `<result>/qc_evaluation.json` и
`<result>/neoqc_qc_report.html`. Отчёт описан в [`html-report.md`](html-report.md).

Русский SVG создаётся из той же фигуры сразу после сохранения английских
файлов: все текстовые элементы переводятся через `scripts/neoqc_i18n.py`,
раскладка пересчитывается (`tight_layout`). При `--formats png` русские SVG не
создаются.

## Манифест

`plots_manifest.json` содержит запись для каждой поддерживаемой пары
«метрика × рид», включая построенные, пропущенные и завершившиеся ошибкой.
Пути в манифесте относительные — от каталога графиков.

```json
{
  "schema_version": 1,
  "theme": "neo-report",
  "formats": ["svg", "png"],
  "locales": ["en", "ru"],
  "figure": {"width_inches": 8.0, "height_inches": 4.5, "aspect_ratio": "16:9",
             "png_width_px": 2400, "png_height_px": 1350, "png_dpi": 300},
  "summary": {"generated": 16, "errors": 0, "skipped": 0},
  "plots": [
    {
      "id": "per_base_quality",
      "read": "R1",
      "title": "Per base sequence quality",
      "source": "per_cycle_R1.tsv",
      "status": "generated",
      "svg": "per_base_quality_R1.svg",
      "png": "per_base_quality_R1.png",
      "localized": {"ru": {"svg": "per_base_quality_R1.ru.svg"}},
      "alt_text": "Per base sequence quality across 38 R1 read groups."
    }
  ]
}
```

| Поле записи | Значение |
|---|---|
| `id`, `read` | идентификатор метрики и рид (`R1`, `R2`) |
| `title` | английское название модуля; перевод берёт отчёт |
| `status` | `generated`, `skipped` или `error` |
| `reason` | причина пропуска или текст ошибки: `source_not_found`, `adapter_analysis_disabled`, сообщение исключения |
| `svg`, `png` | имена файлов для построенных графиков |
| `localized` | необязательно: `{"<язык>": {"svg": "<файл>"}}` для переведённых вариантов |
| `warnings` | необязательно: список проблем с переводом; основной график при этом остаётся `generated` |
| `alt_text` | текстовое описание графика для доступности |

Поля `locales`, `localized` и `warnings` добавлены без смены `schema_version`:
старые потребители их игнорируют, а отчёт без `localized` показывает английские
графики в обоих языках.

Поддерживаемые идентификаторы графиков:

- `per_base_quality`;
- `adapter_content`;
- `per_base_sequence_content`;
- `per_sequence_gc_content`;
- `per_base_n_content`;
- `sequence_length_distribution`;
- `sequence_duplication_levels`;
- `per_sequence_quality`.

Если R2 не передан, его записи получают `source_not_found`. При
`--skip-adapters` записи адаптеров получают `adapter_analysis_disabled`,
остальные графики строятся как обычно. При повторном запуске устаревшие файлы
графиков, включая `*.ru.svg`, удаляются до построения новых.

## Входные TSV

`per_cycle_R1.tsv` / `per_cycle_R2.tsv` содержат `cycle`, `mean_quality`,
`lower_quartile` и `median`. Последние две колонки нужны для
FastQC-совместимой оценки качества по позициям. Старые файлы только со средним
значением строятся, но получают `NOT EVALUATED`, а не выведенный PASS.

`per_sequence_quality_R1.tsv` / `R2.tsv` содержат `mean_quality`, `read_count`
и `read_count_truncate`. Округлённое распределение — собственное представление
NeoQC; усечённое совпадает с разбиением FastQC и используется
FastQC-совместимым профилем. Старые файлы без `read_count_truncate` читаются,
вместо него берётся `read_count`.

`sequence_duplication_levels_R1.tsv` / `R2.tsv` содержат `duplication_level`,
`total_sequences_percent` и `deduplicated_sequences_percent`. Уровни — метки
(`1`, `2`, `>10` …), обе доли — значения от 0 до 100.

Нативный анализ FASTQ также пишет `sequence_duplication_summary_R1.tsv` / `R2.tsv`
— однострочную запись происхождения: идентификатор алгоритма, имя исходного
FASTQ, длину префикса, общее и уникальное число последовательностей и точное
значение `deduplicated_remaining_percent`, которое использует движок статусов.
Доля уровня 1 оставлена только как запасной вариант для ранее импортированных
двухсерийных TSV из FastQC.

Результаты дупликации публикуются как транзакция: в начале запуска NeoQC
удаляет устаревшие артефакты, публикует каждый файл атомарно, сводку — последней,
и только после успеха удаляет маркер `sequence_duplication_R1.incomplete` /
`R2.incomplete`. Пока маркер существует, оценщик отклоняет набор, поэтому
аварийно прерванный запуск не может незаметно подставить неполный профиль.
Метод описан в [`sequence-duplication.md`](sequence-duplication.md).

## Отдельный запуск

```bash
python3 scripts/plot_results.py results/sample01 results/sample01/plots
```

Только один формат и строгий код возврата при ошибках построения:

```bash
python3 scripts/plot_results.py results/sample01 results/sample01/plots \
  --formats svg --strict
```

Явный набор правил:

```bash
python3 scripts/plot_results.py results/sample01 results/sample01/plots \
  --ruleset config/qc_rules/fastqc-compatible-v1.json
```

Оценка без построения графиков:

```bash
python3 scripts/evaluate_qc.py results/sample01
```

Пересборка только HTML-отчёта:

```bash
python3 scripts/generate_qc_report.py results/sample01
python3 scripts/generate_qc_report.py results/sample01 \
  --output results/sample01/sample01_qc.html
```

Ошибка построения графиков не делает недействительными TSV и сводки: C++
считает её предупреждением. Состояние артефакта (`generated`, `skipped`,
`error`) хранится отдельно от QC-статуса и никогда не превращается в QC `FAIL`
(см. [`qc-status-engine.md`](qc-status-engine.md)).

## Визуальный язык

Все графики используют `scripts/plot_style.py`; отдельные функции построения не
задают собственную тему.

- **Шрифты.** IBM Plex Sans для текста, IBM Plex Mono для делений осей. Файлы
  лежат в `assets/fonts/ibm-plex/` и регистрируются в Matplotlib при импорте
  `plot_style`. Если их нет, используется DejaVu Sans / DejaVu Sans Mono.
- **Текст в SVG остаётся текстом** (`svg.fonttype = "none"`) с запасными
  семействами (`'IBM Plex Sans', 'Segoe UI', Arial, sans-serif`). В отчёте текст
  рисуется встроенными веб-шрифтами, поэтому выглядит одинаково в любой системе.
- **Без заголовка внутри фигуры.** Название модуля и рид показывает панель
  отчёта; в метаданных файлов (`Title`) они сохраняются.
- **Размер.** Фигура 8 × 4,5 дюйма (16:9), базовый кегль 11 pt, деления 10 pt.
- **Цвета — только из палитры.** Все цвета берутся из констант `plot_style.py`,
  которые определены токенами `scripts/neoqc_theme.py`. Отчёт заменяет эти цвета
  CSS-переменными, поэтому график перекрашивается в тёмной теме. Цвет вне
  палитры (например, `'red'`) останется неизменным; тест
  `test_russian_chart_variants_and_theme_palette` это запрещает.
- **Нуклеотиды** — как в геномных браузерах: A — зелёный, C — синий, G —
  янтарный, T — красный, N — серый.
- **Строки** на графиках пишутся по-английски; русский перевод — в `CHART` и
  `CHART_PATTERNS` модуля `scripts/neoqc_i18n.py`.
