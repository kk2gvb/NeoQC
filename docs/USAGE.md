# Использование NeoQC

## 1. Общая схема

NeoQC принимает один FASTQ-файл для single-end анализа или пару R1/R2 для paired-end анализа.

```text
FASTQ
  │
  ▼
Проверка структуры
  │
  ▼
Расчёт метрик
  │
  ├── качество
  ├── GC
  ├── N
  ├── длины
  ├── дупликация
  ├── адаптеры
  └── другие метрики
  │
  ▼
TSV-артефакты
  │
  ▼
QC evaluation
  │
  ▼
Графики + HTML report
```

## 2. Single-end

Минимальный запуск:

```bash
./build/neoqc     --r1 sample.fastq.gz     --sample-id sample01     --out results/sample01
```

С графиками и HTML-отчётом:

```bash
./build/neoqc     --r1 sample.fastq.gz     --sample-id sample01     --out results/sample01     --plot
```

## 3. Paired-end

```bash
./build/neoqc     --r1 sample_R1.fastq.gz     --r2 sample_R2.fastq.gz     --sample-id sample01     --out results/sample01     --plot
```

NeoQC проверяет соответствие парных reads и формирует результаты для R1 и R2.

## 4. Параметры командной строки

### `--r1`

Первый FASTQ-файл:

```bash
--r1 sample_R1.fastq.gz
```

### `--r2`

Второй FASTQ-файл для paired-end:

```bash
--r2 sample_R2.fastq.gz
```

### `--sample-id`

Идентификатор образца:

```bash
--sample-id patient01_tumor
```

### `--out`

Каталог результатов:

```bash
--out results/patient01_tumor
```

### `--plot`

Включает построение графиков, QC evaluation и HTML-отчёта:

```bash
--plot
```

### `--skip-adapters`

Отключает расчёт содержания адаптеров:

```bash
--skip-adapters
```

### `--timing`

Выводит время выполнения этапов:

```bash
--timing
```

Пример:

```bash
./build/neoqc     --r1 sample_R1.fastq.gz     --r2 sample_R2.fastq.gz     --sample-id sample01     --out results/sample01     --plot     --timing
```

### `--samples`

Использует CSV-таблицу образцов для проверки или batch-анализа:

```bash
--samples examples/samples.csv
```

### `--help`

```bash
./build/neoqc --help
```

## 5. Batch-анализ

Формат таблицы:

```text
patient_id,sample_id,sample_role,material,r1,r2,platform,library_type,reference
```

Проверить таблицу:

```bash
./build/neoqc --samples examples/samples.csv
```

Запустить batch QC:

```bash
./build/neoqc     --samples examples/samples.csv     --out results     --plot
```

Результаты располагаются по `patient_id` и `sample_id`.

## 6. Выходные данные

Типичный каталог:

```text
results/sample01/
├── *_summary.txt
├── per_cycle_R*.tsv
├── per_sequence_quality_R*.tsv
├── per_base_sequence_content_R*.tsv
├── per_sequence_gc_content_R*.tsv
├── per_base_n_content_R*.tsv
├── sequence_length_distribution_R*.tsv
├── sequence_duplication_levels_R*.tsv
├── sequence_duplication_summary_R*.tsv
├── overrepresented_sequences_R*.tsv
├── adapter_content_R*.tsv
├── qc_evaluation.json
├── neoqc_qc_report.html
└── plots/
```

Для single-end R2-файлы отсутствуют.

## 7. Интерпретация результатов

Подробное описание метрик и графиков:

- [`QC_METRICS.md`](QC_METRICS.md)
- [`FastQ_origin_paired_end_and_plot_legends.md`](FastQ_origin_paired_end_and_plot_legends.md)

HTML-отчёт предназначен для просмотра человеком, а `qc_evaluation.json` — для машинной обработки и проверки статусов.

## 8. QC-статусы

```text
рассчитанная метрика
        │
        ▼
наблюдение
        │
        ▼
versioned ruleset
        │
        ▼
PASS / WARNING / FAIL
```

Статус определяется правилами, а не самим графиком.

## 9. Повторный запуск

NeoQC изолирует результаты отдельного запуска во временном каталоге и публикует новый набор результатов после успешного завершения анализа.

Это предотвращает смешивание артефактов разных запусков, например:

```text
paired-end → single-end
```

или:

```text
с поиском адаптеров → --skip-adapters
```

При повторном запуске можно использовать тот же `--out`, чтобы заменить предыдущий результат.

## 10. Дупликация

NeoQC рассчитывает дупликацию по первым 50 нуклеотидам каждого read.

Каждый уникальный ключ учитывается на протяжении всего входного файла.

Подробности:

[`sequence-duplication.md`](sequence-duplication.md)

## 11. Производительность

Для измерения времени:

```bash
./build/neoqc     --r1 sample_R1.fastq.gz     --r2 sample_R2.fastq.gz     --sample-id sample01     --out results/sample01     --timing
```

NeoQC использует OpenMP для параллельной обработки FASTQ.

## 12. Типовой рабочий сценарий

### Шаг 1. Получение FASTQ

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

### Шаг 2. Первичный QC

```bash
./build/neoqc     --r1 sample_R1.fastq.gz     --r2 sample_R2.fastq.gz     --sample-id sample01     --out results/sample01     --plot
```

### Шаг 3. Просмотр отчёта

```text
results/sample01/neoqc_qc_report.html
```

### Шаг 4. Предварительная обработка

При необходимости выполняется отдельным инструментом trimming/filtering.

NeoQC эти операции не выполняет.

### Шаг 5. Повторный QC

После preprocessing можно снова запустить NeoQC на полученных FASTQ.

### Шаг 6. Передача данных дальше

После контроля качества данные могут использоваться на последующих этапах NGS-пайплайна.

## 13. Ошибки входных данных

NeoQC проверяет, среди прочего:

- структуру FASTQ;
- начало записи;
- наличие `+`;
- соответствие длины sequence и quality;
- допустимые символы качества;
- допустимые основания;
- целостность gzip;
- соответствие paired-end reads.

При ошибке входные данные не следует интерпретировать как прошедшие QC.

## 14. Что NeoQC не делает

NeoQC не является полным NGS-пайплайном.

Он не выполняет:

```text
trimming
filtering
alignment
variant calling
expression quantification
HLA typing
neoepitope prediction
```

Его задача — контроль качества FASTQ и формирование проверяемого технического отчёта.
