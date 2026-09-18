# NeoQC — окончательный план CUDA-миграции и оптимизации

## 1. Цель

Добавить в NeoQC второй вычислительный backend — CUDA — без удаления существующей реализации OpenMP.

Итоговая модель:

```text
                    NeoQC
                      │
                --compute
                 /        \
              omp          cuda
               │             │
          CPU / OMP      OMP + CUDA
               │             │
               └──────┬──────┘
                      │
                common results
                      │
             QC / TSV / JSON / plots
```

### Режимы

#### `--compute omp`

Существующая реализация NeoQC:

```text
FASTQ
  ↓
OpenMP / CPU
  ↓
merge
  ↓
results
```

Этот режим сохраняется как основной fallback и reference implementation.

#### `--compute cuda`

Гибридный режим:

```text
FASTQ
  ↓
OpenMP / CPU
  ├── reading
  ├── decompression
  ├── batch preparation
  └── CPU-only operations
          ↓
        CUDA
  ├── тяжелые массово-параллельные расчеты
  └── GPU reductions
          ↓
        OMP / CPU
  ├── final aggregation
  ├── QC evaluation
  └── output
```

Важно: выбор `cuda` **не означает перенос всего проекта на GPU**. В CUDA-режиме OpenMP продолжает использоваться для CPU-частей pipeline.

---

# 2. Ограничения проекта

## 2.1. K-mer

**K-mer counting в NeoQC не входит.**

Ранее было принято решение не добавлять эту метрику. Поэтому:

- K-mer не является функциональностью NeoQC;
- K-mer не является CUDA-задачей;
- K-mer отсутствует в таблицах миграции;
- K-mer не должен появляться в итоговом pipeline.

## 2.2. Существующие CLI-режимы

Существующие параметры NeoQC не изменяются:

```text
--r1
--r2
--sample-id
--out
--plot
--skip-adapters
--timing
--samples
```

Добавляется только выбор вычислительного backend:

```text
--compute omp
--compute cuda
```

Точное имя флага можно окончательно определить при реализации.

---

# 3. Полная карта компонентов NeoQC

| Компонент / операция | Сейчас | `omp` | `cuda` | Стратегия |
|---|---|---:|---:|---|
| CLI | CPU | OMP/CPU | OMP/CPU | Оставить |
| Sample sheet parsing | CPU | OMP/CPU | OMP/CPU | Оставить |
| Sample sheet validation | CPU | OMP/CPU | OMP/CPU | Оставить |
| FASTQ reading | CPU | OMP/CPU | OMP/CPU | Оставить |
| gzip/decompression | CPU | OMP/CPU | OMP/CPU | Оставить на первом этапе |
| Batch preparation | CPU | OMP/CPU | OMP/CPU | Оставить |
| Read representation | C++ structures | OMP/CPU | Host + Device buffers | Переделать для CUDA |
| Per-base quality | OMP | OMP | CUDA | **P0: переносить** |
| Per-cycle quality | OMP | OMP | CUDA | **P0: переносить** |
| N content | OMP | OMP | CUDA | **P0: переносить** |
| GC statistics | OMP | OMP | CUDA | **P0: переносить** |
| Per-sequence quality | OMP | OMP | CUDA | **P1: переносить** |
| Per-sequence GC | OMP | OMP | CUDA | **P1: переносить** |
| Adapter statistics | OMP | OMP | CUDA/OMP | **P2: точечный перенос** |
| Sequence duplication | OMP | OMP | CUDA/OMP | **P2: отдельное проектирование** |
| Local statistics | OMP thread-local | OMP | GPU counters/reduction | Переделать |
| Merge | OMP/CPU | OMP/CPU | OMP/CPU + GPU reduction | **Отдельно оптимизировать** |
| Final statistics | CPU | OMP/CPU | OMP/CPU | Оставить |
| QC evaluation | Python/CPU | CPU | CPU | Оставить |
| TSV generation | CPU | CPU | CPU | Оставить |
| JSON generation | CPU | CPU | CPU | Оставить |
| Plot generation | Python | CPU | CPU | Оставить |
| K-mer counting | Нет | Нет | Нет | **Не добавлять** |

---

# 4. Что точно оставить на OMP/CPU

Следующие операции не являются приоритетом CUDA-миграции.

## 4.1. Ввод-вывод

```text
FASTQ reading
gzip/decompression
file handling
```

Причина: это преимущественно I/O и подготовка данных. Перенос на GPU сам по себе не означает ускорение.

## 4.2. Управление pipeline

```text
CLI
sample sheet
configuration
orchestration
error handling
```

GPU здесь не требуется.

## 4.3. Формирование результатов

```text
TSV
JSON
summary
QC evaluation
plot generation
```

Эти операции работают с относительно небольшим объемом уже агрегированных данных.

## 4.4. Другие небольшие CPU-операции

Любая операция остается на OMP/CPU, если профилирование показывает, что:

```text
GPU kernel + H2D + D2H
```

дороже самой CPU-операции.

---

# 5. Что точечно переделывать на CUDA

## P0 — основной кандидат

### 5.1. Per-base / Per-cycle quality

Причины для CUDA:

- обрабатываются огромные количества bases;
- операции над reads хорошо распараллеливаются;
- результат можно представить компактными счетчиками;
- возможен GPU reduction.

Целевая схема:

```text
reads
  ↓
CUDA threads
  ↓
quality counters
  ↓
reduction
  ↓
compact result
```

### 5.2. N content

Аналогично:

```text
base == 'N'
    ↓
parallel counting
    ↓
per-cycle counters
```

### 5.3. GC statistics

```text
G/C detection
    ↓
parallel counting
    ↓
per-cycle / per-sequence aggregation
```

---

# 6. P1 — следующие CUDA-кандидаты

## 6.1. Per-sequence quality

Перенос после стабилизации P0.

Особое внимание:

- histogram;
- размеры структур;
- reduction;
- различия результатов OMP/CUDA.

## 6.2. Per-sequence GC

Переносить после общей инфраструктуры batch и reduction.

---

# 7. P2 — сложные операции

## 7.1. Adapter statistics

Возможна CUDA-реализация, но сначала определить:

- какие операции выполняются для каждого read;
- какие структуры используются;
- есть ли branching;
- как выполняется aggregation;
- насколько велика стоимость передачи данных.

Цель — не переносить весь adapter pipeline автоматически, а перенести именно тяжелую часть.

## 7.2. Sequence duplication

Это отдельная задача.

Причина:

```text
duplication
    ↓
большое количество последовательностей
    ↓
hash/map structures
    ↓
aggregation / merge
```

GPU-реализация может потребовать принципиально другой структуры данных.

Поэтому сначала:

1. профилирование;
2. анализ памяти;
3. анализ текущего merge;
4. оценка GPU-варианта;
5. только затем реализация CUDA.

---

# 8. Проблема `Killed`

Это отдельный P0.

Наблюдаемая ситуация:

```text
Processed 100%
      ↓
работа с reads завершена
      ↓
merge / aggregation / finalization
      ↓
Killed
```

Пока нельзя утверждать, что причиной является именно `merge()`.

Возможные точки:

```text
last batch
    ↓
temporary allocations
    ↓
merge
    ↓
duplication merge
    ↓
sorting
    ↓
final statistics
    ↓
output preparation
```

## 8.1. Сначала подтвердить OOM

Запуск:

```bash
/usr/bin/time -v ./build/neoqc ...
```

Проверить:

```text
Maximum resident set size
```

После падения:

```bash
dmesg -T | grep -i -E 'oom|killed process|out of memory'
```

или:

```bash
journalctl -k | grep -i -E 'oom|killed process'
```

## 8.2. Что искать в коде

Проверить:

- хранение всех reads;
- большие `std::vector`;
- `unordered_map` / `map`;
- thread-local analyzers;
- временные копии;
- `merge()` с созданием второго большого объекта;
- duplication structures;
- sorting;
- преобразование результатов;
- одновременное существование local + global + temporary structures.

## 8.3. Цель оптимизации merge

Нежелательно:

```text
local results
     ↓
создание огромного temporary result
     ↓
copy
     ↓
merge
```

Предпочтительно:

```text
local result
     ↓
incremental merge
     ↓
global result
     ↓
release local result
```

Если возможно, merge должен быть потоковым/инкрементальным.

---

# 9. Архитектура данных для CUDA

Текущие C++-структуры, содержащие `std::string` и другие сложные объекты, не следует напрямую передавать GPU.

Нужен общий batch:

```cpp
struct ReadBatch {
    std::vector<char> sequences;
    std::vector<char> qualities;
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> lengths;
};
```

Для CUDA:

```text
ReadBatch
   ↓
Host buffers
   ↓
Device buffers
```

Например:

```text
sequence_data
quality_data
offsets
lengths
```

Это позволяет использовать один логический batch для обоих режимов.

---

# 10. Общая backend-архитектура

Целевая модель:

```text
                     NeoQC
                       │
                 --compute
                  /       \
               omp         cuda
                │            │
                ▼            ▼
             OMP path    hybrid path
                │            │
                │       ┌────┴─────┐
                │       │          │
                │      OMP       CUDA
                │       │          │
                └───────┴──────────┘
                        │
                  common results
                        │
                QC / reports
```

Лучше выделить общий интерфейс вычислительного backend.

Например концептуально:

```cpp
class IAnalyzerBackend {
public:
    virtual void processBatch(const ReadBatch& batch) = 0;
    virtual QualityStats getResults() = 0;
    virtual ~IAnalyzerBackend() = default;
};
```

Реализации:

```text
OMPBackend
CudaBackend
```

При этом `CudaBackend` внутри может использовать OMP для CPU-частей.

---

# 11. CUDA kernel strategy

Не переносить OpenMP-код буквально.

Текущая модель:

```text
OMP thread
    ↓
local QualityAnalyzer
    ↓
processRecord()
    ↓
merge()
```

CUDA-модель:

```text
GPU threads
    ↓
compact counters
    ↓
parallel reduction
    ↓
small result
    ↓
CPU/OMP
```

Цель:

> как можно больше вычислений выполнять на GPU и как можно меньше промежуточных данных возвращать на CPU.

---

# 12. Reduction

Reduction является одной из ключевых задач CUDA-миграции.

Нежелательно:

```text
каждый GPU thread
    ↓
большой analyzer
    ↓
Device → Host
    ↓
огромный набор результатов
    ↓
CPU merge
```

Предпочтительно:

```text
GPU threads
     ↓
local counters
     ↓
block reduction
     ↓
global reduction
     ↓
compact result
     ↓
CPU
```

Это одновременно:

- снижает объем D2H;
- снижает число объектов;
- уменьшает CPU merge;
- потенциально снижает пиковое потребление RAM.

---

# 13. Batch size

Необходимо отдельно протестировать размеры batch:

```text
10k
50k
100k
250k
500k
1M reads
```

Для каждого измерять:

| Показатель | Измерять |
|---|---:|
| Total time | Да |
| OMP time | Да |
| CUDA kernel time | Да |
| H2D | Да |
| D2H | Да |
| Merge | Да |
| RAM | Да |
| VRAM | Да |
| GPU utilization | Да |

Оптимальный batch выбирается по совокупности производительности и памяти.

---

# 14. Передача данных

Измерять отдельно:

```text
T_total =
    T_read
  + T_prepare
  + T_H2D
  + T_kernel
  + T_D2H
  + T_merge
  + T_output
```

Нельзя оценивать CUDA только по времени kernel.

Например:

```text
kernel = 1 s
H2D + D2H = 8 s
```

не означает ускорение всего NeoQC.

---

# 15. Пошаговая миграция

## P0 — диагностика текущей реализации

- [ ] Запустить `/usr/bin/time -v`.
- [ ] Проверить kernel log на OOM.
- [ ] Найти точную фазу `Killed`.
- [ ] Замерить RAM перед merge.
- [ ] Замерить RAM во время merge.
- [ ] Замерить RAM после merge.
- [ ] Исследовать duplication merge.
- [ ] Найти временные копии.
- [ ] Найти структуры, удерживающие данные после `Processed 100%`.

## P1 — backend infrastructure

- [ ] Определить `--compute omp`.
- [ ] Определить `--compute cuda`.
- [ ] Оставить OMP backend без изменения поведения.
- [ ] Выделить общий backend interface.
- [ ] Выделить общий `ReadBatch`.
- [ ] Подготовить CMake для optional CUDA.
- [ ] CUDA не должна быть обязательной для OMP-сборки.

## P2 — CUDA foundation

- [ ] Добавить `.cu` / `.cuh`.
- [ ] Инициализация CUDA.
- [ ] Проверка устройства.
- [ ] Host → Device.
- [ ] Device → Host.
- [ ] Первый CUDA kernel.
- [ ] Первый GPU reduction.
- [ ] Сравнение с OMP.

## P3 — базовые метрики

- [ ] Per-base quality.
- [ ] Per-cycle quality.
- [ ] N content.
- [ ] GC statistics.

После каждого пункта:

```text
OMP result
    vs
CUDA result
```

## P4 — sequence-level metrics

- [ ] Per-sequence quality.
- [ ] Per-sequence GC.
- [ ] Проверка histogram/reduction.
- [ ] Проверка больших read counts.
- [ ] Проверка mixed-length reads.

## P5 — сложные метрики

- [ ] Adapter statistics — точечный перенос.
- [ ] Sequence duplication — отдельное проектирование.
- [ ] Оценка необходимости GPU merge.

## P6 — memory optimization

- [ ] Уменьшить local result structures.
- [ ] Уменьшить temporary allocations.
- [ ] Сделать incremental merge там, где возможно.
- [ ] Перенести aggregation на GPU там, где это уменьшает объем результатов.
- [ ] Проверить peak RAM.
- [ ] Проверить peak VRAM.

## P7 — benchmark

- [ ] OMP benchmark.
- [ ] CUDA benchmark.
- [ ] Hybrid benchmark.
- [ ] Разные batch sizes.
- [ ] Разные размеры FASTQ.
- [ ] Слабый компьютер.
- [ ] Мощный компьютер с NVIDIA GPU.

## P8 — стабилизация

- [ ] Все существующие тесты.
- [ ] CUDA-specific tests.
- [ ] Cross-backend regression tests.
- [ ] Малые FASTQ.
- [ ] Большие FASTQ.
- [ ] Single-end.
- [ ] Paired-end.
- [ ] Mixed-length reads.
- [ ] Проверка `--skip-adapters`.
- [ ] Проверка `--timing`.
- [ ] Проверка `--samples`.
- [ ] Проверка запуска без CUDA через `--compute omp`.

---

# 16. Критерий завершения

CUDA-миграция считается технически завершенной, когда:

1. `--compute omp` продолжает работать.
2. `--compute cuda` запускает гибридный OMP+CUDA pipeline.
3. Существующие CLI-режимы не меняют поведения.
4. Основные массово-параллельные QC-метрики работают на CUDA.
5. Результаты CUDA совпадают с OMP в пределах допустимых численных различий.
6. CUDA не требует переноса I/O и отчетности на GPU.
7. Пиковое потребление памяти контролируется.
8. Причина `Killed` установлена и устранена либо точно локализована.
9. CUDA benchmark измеряет H2D, kernel, D2H и merge отдельно.
10. Необходимость дальнейшего переноса Adapter/Duplication определяется профилированием, а не предположением.
11. K-mer counting в проект не добавляется.
12. На машине без NVIDIA GPU остается полностью рабочий `--compute omp`.

---

# 17. Итоговая целевая схема

```text
                         NeoQC
                           │
                     --compute
                     /          \
                  omp            cuda
                   │               │
                   ▼               ▼
             ┌──────────┐    ┌──────────────┐
             │ OMP/CPU  │    │ OMP/CPU      │
             │          │    │ FASTQ        │
             │ all QC   │    │ preparation  │
             │ metrics  │    │ CPU-only     │
             └────┬─────┘    └──────┬───────┘
                  │                 │
                  │                 ▼
                  │           ┌──────────────┐
                  │           │    CUDA      │
                  │           │              │
                  │           │ quality      │
                  │           │ GC           │
                  │           │ N            │
                  │           │ sequence QC  │
                  │           │ reductions   │
                  │           └──────┬───────┘
                  │                  │
                  └─────────┬────────┘
                            ▼
                     common results
                            │
                            ▼
                    merge / QC evaluation
                            │
                     TSV / JSON / plots
```

Главный принцип проекта:

> **Не «переписать NeoQC на CUDA», а добавить CUDA как второй backend и постепенно перенести на GPU те вычислительные участки, где это дает измеримый выигрыш. OpenMP при этом остается полноценным CPU backend и используется также внутри гибридного CUDA-режима.**
