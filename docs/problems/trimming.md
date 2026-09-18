# Trimming в NeoQC

## Цель

Добавить в NeoQC preprocessing FASTQ-файлов с возможностью удаления адаптеров,
низкокачественных участков и слишком коротких reads.

Обработка должна сохранять исходные FASTQ-файлы и создавать очищенные FASTQ
как отдельный результат.

Общий pipeline:

```text
FASTQ
  ↓
QC before trimming
  ↓
Trimming
  ├── fixed trimming
  ├── quality trimming
  ├── adapter trimming
  ├── polyG/polyX trimming
  └── length filtering
  ↓
QC after trimming
  ↓
clean FASTQ
  ↓
alignment / дальнейший pipeline
```

---

# Design decisions

## 1. Отдельный модуль trimming

Trimming не должен засорять основные директории `src/` и `include/`.

Вся функциональность должна находиться в отдельном модуле:

```text
include/trimming/
src/trimming/
tests/trimming/
```

Предполагаемая структура:

```text
include/trimming/
├── trim_config.h
├── trim_result.h
├── trim_stats.h
├── quality_trimmer.h
├── adapter_trimmer.h
└── trimmer.h

src/trimming/
├── quality_trimmer.cpp
├── adapter_trimmer.cpp
└── trimmer.cpp

tests/trimming/
└── ...
```

### Правила

- публичный API trimming находится в `include/trimming/`;
- реализация находится в `src/trimming/`;
- тесты находятся в `tests/trimming/`;
- не добавлять файлы trimming непосредственно в `src/` или `include/`;
- не смешивать алгоритмы trimming с существующей логикой QC.

---

## 2. `--trim` — основной feature flag

Флаг должен появиться в самом начале разработки функции, сразу после
создания структуры модуля.

```bash
neoqc --r1 input.fastq.gz
```

должен работать так же, как до добавления trimming.

```bash
neoqc --r1 input.fastq.gz --trim
```

включает trimming pipeline.

### Принцип

```text
--trim
   │
   ▼
trimming/
   │
   ├── fixed trimming
   ├── quality trimming
   ├── adapter trimming
   └── length filtering
```

`--trim` является главным переключателем функции.

Остальные CLI-параметры управляют отдельными этапами trimming.

---

## 3. Backward compatibility

Добавление trimming не должно менять существующее поведение NeoQC,
если `--trim` не указан.

Обязательный regression test:

```text
NeoQC без --trim
        ↓
старое поведение
        ↓
старые тесты должны проходить
```

---

## 4. Разделение trimming и filtering

Trimming и filtering являются разными операциями.

```text
TRIMMING
↓
изменение read

FILTERING
↓
решение оставить / удалить read
```

Например:

```text
150 bp
 ↓
adapter trimming
 ↓
120 bp
 ↓
quality trimming
 ↓
100 bp
 ↓
length filtering
 ↓
PASS
```

или:

```text
150 bp
 ↓
trimming
 ↓
35 bp
 ↓
min_length = 50
 ↓
DISCARD
```

Причина удаления read не должна смешиваться с причиной trimming.

---

# P0 — Infrastructure

## TRIM-001 — Создать модуль trimming

Создать:

```text
include/trimming/
src/trimming/
tests/trimming/
```

### Требования

- изолировать trimming от остальных модулей NeoQC;
- подготовить подключение нового модуля в CMake;
- существующая сборка NeoQC не должна ломаться;
- существующие тесты должны продолжать проходить.

---

## TRIM-002 — Добавить `--trim`

Добавить основной feature flag:

```text
--trim
```

### Поведение

Без флага:

```bash
neoqc --r1 input.fastq.gz
```

работает без trimming.

С флагом:

```bash
neoqc --r1 input.fastq.gz --trim
```

запускается trimming pipeline.

### Требования

- добавить параметр в CLI;
- передавать состояние флага в конфигурацию;
- пока конкретные алгоритмы не реализованы, флаг не должен ломать pipeline;
- сохранить backward compatibility.

### Тесты

Проверить оба режима:

```text
--trim отсутствует
--trim присутствует
```

---

## TRIM-003 — Создать TrimConfig

Файл:

```text
include/trimming/trim_config.h
```

Содержит настройки trimming.

### Предварительные параметры

```text
enabled

trim_front
trim_tail

cut_front
cut_tail
cut_right

quality_threshold
window_size

adapter_trimming
adapter_sequence

min_length
```

Для paired-end предусмотреть отдельные параметры R1/R2 там,
где это необходимо.

---

## TRIM-004 — Создать TrimResult

Файл:

```text
include/trimming/trim_result.h
```

Описывает результат обработки одного read.

### Минимальные поля

```text
original_length
final_length

trimmed_front
trimmed_tail

adapter_found
adapter_position

passed
discard_reason
```

---

## TRIM-005 — Создать TrimStats

Файл:

```text
include/trimming/trim_stats.h
```

Содержит агрегированную статистику trimming.

### Минимальные поля

```text
total_reads
passed_reads
discarded_reads

bases_before
bases_after
bases_trimmed

adapter_trimmed_reads
quality_trimmed_reads
too_short_reads
```

---

## TRIM-006 — Создать QualityTrimmer

Файлы:

```text
include/trimming/quality_trimmer.h
src/trimming/quality_trimmer.cpp
```

Ответственность:

- fixed front trimming;
- fixed tail trimming;
- quality front trimming;
- quality tail trimming;
- sliding-window trimming.

Алгоритмы не должны быть реализованы в `Trimmer`.

---

## TRIM-007 — Создать AdapterTrimmer

Файлы:

```text
include/trimming/adapter_trimmer.h
src/trimming/adapter_trimmer.cpp
```

Ответственность:

- adapter sequence;
- adapter matching;
- mismatches;
- PE overlap trimming.

---

## TRIM-008 — Создать Trimmer

Файлы:

```text
include/trimming/trimmer.h
src/trimming/trimmer.cpp
```

`Trimmer` является верхнеуровневым обработчиком.

Он объединяет:

```text
QualityTrimmer
AdapterTrimmer
Length filtering
TrimStats
```

`Trimmer` отвечает за порядок выполнения операций, но не содержит
реализацию отдельных алгоритмов.

---

# P0 — Algorithms

## TRIM-009 — Fixed front trimming

Добавить удаление фиксированного количества bases с 5'-конца read.

Пример:

```text
--trim-front 10
```

### Тесты

- trimming 0 bases;
- trimming 1 base;
- trimming N bases;
- trimming больше длины read;
- paired-end R1/R2.

---

## TRIM-010 — Fixed tail trimming

Добавить удаление фиксированного количества bases с 3'-конца read.

Пример:

```text
--trim-tail 10
```

### Тесты

- trimming 0 bases;
- trimming 1 base;
- trimming N bases;
- trimming больше длины read;
- paired-end.

---

## TRIM-011 — Quality front trimming

Удалять низкокачественный участок с начала read.

```text
5' → 3'

[LOW][LOW][LOW][GOOD][GOOD][GOOD]
  ↓
 удалить
             ↓
          сохранить
```

### Параметры

```text
quality_threshold
window_size
```

### Тесты

- полностью качественный read;
- низкое качество в начале;
- низкое качество в середине;
- низкое качество в конце;
- несколько низкокачественных участков.

---

## TRIM-012 — Quality tail trimming

Удалять низкокачественный участок с конца read.

```text
[GOOD][GOOD][GOOD][LOW][LOW][LOW]
                     ↓
                   удалить
```

### Тесты

- плохой tail;
- хороший tail;
- read полностью низкого качества;
- короткий read после trimming.

---

## TRIM-013 — Sliding-window trimming

Добавить trimming по скользящему окну.

### Параметры

```text
window_size
quality_threshold
```

### Требования

Использовать rolling sum / rolling average, чтобы не пересчитывать
среднее качество окна с нуля для каждой позиции.

### Тесты

- окно полностью хорошего качества;
- окно с низким качеством;
- низкое качество в середине;
- низкое качество в конце;
- `window_size = 1`;
- `window_size = read_length`;
- `window_size > read_length`.

---

## TRIM-014 — Adapter sequence matching

Добавить поиск заданной adapter sequence в read.

Пример:

```text
--adapter-sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
```

После обнаружения adapter:

```text
read
───────────────────────────────
              adapter
                ↓
────────────────┬──────────────
                ↓
              trim
```

### Тесты

Создать synthetic reads:

```text
NO_ADAPTER
HALF_ADAPTER
FULL_ADAPTER
```

Проверить:

- adapter не найден;
- adapter найден частично;
- adapter найден полностью.

---

## TRIM-015 — Adapter mismatch handling

Разрешить небольшое количество несовпадений между read и adapter.

### Необходимо определить

- максимальное количество mismatches;
- минимальную длину совпадения;
- критерий принятия adapter match.

### Тесты

- exact match;
- 1 mismatch;
- несколько mismatches;
- слишком короткое совпадение;
- случайное совпадение.

---

## TRIM-016 — PE overlap adapter trimming

Для paired-end reads реализовать поиск adapter через overlap R1/R2.

```text
R1 ────────────────────>
       <──────────────── R2

          overlap
             ↓
       определить insert
             ↓
      удалить adapter
```

### Требования

Для пары reads:

1. определить overlap;
2. проверить overlap;
3. определить границу insert;
4. обрезать adapter;
5. сохранить синхронизацию R1/R2.

### Тесты

- normal PE;
- короткий insert;
- полный overlap;
- отсутствие overlap;
- adapter только в R1;
- adapter только в R2;
- adapter в обоих reads.

---

## TRIM-017 — Minimum length filtering

После trimming удалять reads, которые стали слишком короткими.

Параметр:

```text
--min-length N
```

Пример:

```text
original = 150 bp
trimmed  = 42 bp
min      = 50 bp

→ read discarded
```

### Тесты

- длина > min;
- длина = min;
- длина < min;
- read стал короче после trimming.

---

# P1 — Statistics

## TRIM-018 — Collect trimming statistics

Собирать:

```text
total_reads
passed_reads
discarded_reads

bases_before
bases_after
bases_trimmed

adapter_trimmed_reads
quality_trimmed_reads
polyG_trimmed_reads
polyX_trimmed_reads

too_short_reads
```

---

## TRIM-019 — Adapter position statistics

Сохранять распределение позиции найденного adapter.

Пример:

```text
adapter_position
----------------
101 : 1234
102 : 5321
103 : 8123
104 : 421
```

Эти данные могут использоваться для QC-графика.

---

## TRIM-020 — Generate trimming report

Создавать отдельный trimming report.

Пример:

```json
{
  "reads": {
    "total": 1000000,
    "passed": 982341,
    "discarded": 17659
  },
  "bases": {
    "before": 150000000,
    "after": 141238421,
    "trimmed": 8761579
  }
}
```

Формат отчёта необходимо согласовать с существующей системой
результатов NeoQC.

---

# P1 — CLI parameters

## TRIM-021 — Добавить CLI параметры trimming

Добавить параметры для управления отдельными этапами trimming.

Предварительный набор:

```text
--trim-front
--trim-tail

--cut-front
--cut-tail
--cut-right

--quality-threshold
--window-size

--adapter-sequence
--adapter-fasta

--min-length
```

Названия параметров могут быть изменены во время реализации.

---

# P1 — Output

## TRIM-022 — Добавить output для trimmed FASTQ

Не перезаписывать исходный FASTQ.

Предлагаемая структура:

```text
results/
└── SAMPLE/
    ├── input/
    │   ├── R1.fastq.gz
    │   └── R2.fastq.gz
    │
    ├── trimmed/
    │   ├── R1.trimmed.fastq.gz
    │   └── R2.trimmed.fastq.gz
    │
    └── qc/
```

### Требование

Original FASTQ остаётся неизменным.

---

## TRIM-023 — QC before/after trimming

Добавить два состояния QC:

```text
QC BEFORE
   ↓
TRIMMING
   ↓
QC AFTER
```

Сравнивать:

- mean quality;
- per-cycle quality;
- adapter content;
- GC;
- N content;
- sequence length;
- duplication;
- количество reads.

---

# P1 — Paired-end consistency

## TRIM-024 — Проверить синхронизацию R1/R2

После trimming/filtering R1 и R2 должны оставаться синхронизированными.

До:

```text
R1_001 ↔ R2_001
R1_002 ↔ R2_002
R1_003 ↔ R2_003
```

После удаления пары:

```text
R1_001 ↔ R2_001
R1_003 ↔ R2_003
```

Нельзя получить:

```text
R1_001 ↔ R2_002
```

### Тесты

Проверить filtering/trimming на synthetic PE dataset.

---

# P1 — Extended trimming

## TRIM-025 — PolyG trimming

Добавить удаление polyG tail.

```text
ACGTACGTACGTGGGGGGGGGG
                  ↑
                trim
```

### Параметр

```text
polyG_min_length
```

### Тесты

- короткий polyG;
- длинный polyG;
- polyG в середине;
- polyG в конце.

---

## TRIM-026 — PolyX trimming

Добавить trimming гомополимерных хвостов:

```text
AAAAAA
CCCCCC
GGGGGG
TTTTTT
```

### Параметр

```text
polyX_min_length
```

### Требование

Trim выполняется только для terminal tail.

---

# P2 — Extended PE processing

## TRIM-027 — PE overlap correction

Использовать overlap R1/R2 для коррекции mismatched bases.

---

## TRIM-028 — PE read merging

Добавить возможность объединения overlapping R1/R2.

---

## TRIM-029 — UMI preprocessing

Добавить поддержку UMI preprocessing.

---

## TRIM-030 — Adapter FASTA

Добавить загрузку нескольких adapter sequences из FASTA.

---

# Tests

## TRIM-031 — Synthetic FASTQ dataset

Создать тестовый набор:

```text
clean.fastq
adapter.fastq
low_quality.fastq
polyG.fastq
polyX.fastq
short.fastq

PE_clean_R1.fastq
PE_clean_R2.fastq
PE_adapter_R1.fastq
PE_adapter_R2.fastq
```

---

## TRIM-032 — Unit tests for QualityTrimmer

Покрыть тестами:

- fixed front;
- fixed tail;
- quality front;
- quality tail;
- sliding-window;
- граничные случаи.

---

## TRIM-033 — Unit tests for AdapterTrimmer

Покрыть тестами:

- exact match;
- partial match;
- mismatch;
- отсутствие adapter;
- PE overlap;
- fallback.

---

## TRIM-034 — Unit tests for Trimmer

Проверить правильный порядок операций:

```text
fixed trim
    ↓
quality trim
    ↓
adapter trim
    ↓
length filtering
```

---

## TRIM-035 — Boundary tests

Обязательно проверить:

```text
read length = 0
read length = 1

trim = 0
trim = read_length
trim > read_length

quality threshold = 0
quality threshold = maximum

window = 1
window = read_length
window > read_length

min_length = 0
min_length = read_length
min_length > read_length
```

---

## TRIM-036 — Paired-end tests

Проверить:

- одинаковое количество R1/R2;
- сохранение read IDs;
- удаление пары целиком при filtering;
- trimming R1/R2;
- adapter trimming;
- overlap analysis.

---

## TRIM-037 — CLI tests

Проверить:

```text
без --trim
с --trim
```

и отдельные CLI-параметры trimming.

Особенно важно проверить, что:

```bash
neoqc --r1 input.fastq.gz
```

не меняет существующее поведение NeoQC.

---

## TRIM-038 — Regression tests

После добавления trimming убедиться, что существующий QC pipeline
не изменился при выключенном `--trim`.

```text
NeoQC без --trim
        ↓
старые результаты
        ↓
старые тесты должны проходить
```

---

# Итоговый pipeline

После выполнения P0:

```text
                    ┌──────────────┐
FASTQ ─────────────►│ QC BEFORE    │
                    └──────┬───────┘
                           ↓
                    ┌──────────────┐
                    │ FIXED TRIM   │
                    └──────┬───────┘
                           ↓
                    ┌──────────────┐
                    │ QUALITY TRIM │
                    └──────┬───────┘
                           ↓
                    ┌──────────────┐
                    │ ADAPTER TRIM │
                    └──────┬───────┘
                           ↓
                    ┌──────────────┐
                    │ LENGTH FILTER│
                    └──────┬───────┘
                           ↓
                    ┌──────────────┐
                    │ QC AFTER     │
                    └──────┬───────┘
                           ↓
                    trimmed FASTQ
                           ↓
                       ALIGNMENT
```

---

# Приоритет

## P0 — Базовая версия

- [ ] TRIM-001 — Создать модуль trimming
- [ ] TRIM-002 — Добавить `--trim`
- [ ] TRIM-003 — Создать TrimConfig
- [ ] TRIM-004 — Создать TrimResult
- [ ] TRIM-005 — Создать TrimStats
- [ ] TRIM-006 — Создать QualityTrimmer
- [ ] TRIM-007 — Создать AdapterTrimmer
- [ ] TRIM-008 — Создать Trimmer
- [ ] TRIM-009 — Fixed front trimming
- [ ] TRIM-010 — Fixed tail trimming
- [ ] TRIM-011 — Quality front trimming
- [ ] TRIM-012 — Quality tail trimming
- [ ] TRIM-013 — Sliding-window trimming
- [ ] TRIM-014 — Adapter sequence matching
- [ ] TRIM-015 — Adapter mismatch handling
- [ ] TRIM-016 — PE overlap adapter trimming
- [ ] TRIM-017 — Minimum length filtering

## P1 — Полноценный preprocessing

- [ ] TRIM-018 — Collect trimming statistics
- [ ] TRIM-019 — Adapter position statistics
- [ ] TRIM-020 — Generate trimming report
- [ ] TRIM-021 — CLI параметры trimming
- [ ] TRIM-022 — Trimmed FASTQ output
- [ ] TRIM-023 — QC before/after
- [ ] TRIM-024 — PE synchronization
- [ ] TRIM-025 — PolyG trimming
- [ ] TRIM-026 — PolyX trimming

## P2 — Расширение

- [ ] TRIM-027 — PE overlap correction
- [ ] TRIM-028 — PE read merging
- [ ] TRIM-029 — UMI preprocessing
- [ ] TRIM-030 — Adapter FASTA

## Tests

- [ ] TRIM-031 — Synthetic FASTQ dataset
- [ ] TRIM-032 — QualityTrimmer tests
- [ ] TRIM-033 — AdapterTrimmer tests
- [ ] TRIM-034 — Trimmer tests
- [ ] TRIM-035 — Boundary tests
- [ ] TRIM-036 — Paired-end tests
- [ ] TRIM-037 — CLI tests
- [ ] TRIM-038 — Regression tests
