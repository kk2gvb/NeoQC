# NeoQC — план оптимизации OpenMP и внедрения CUDA с учётом тримминга

**База анализа:** `kk2gvb/NeoQC`, ветка `optimization_and_cuda`, коммит `ac17564b7996995feee2946e94d30998ec8026f8` (29.10.2026).  
**Статус:** технический план; ускорение и OOM не подтверждены локальными замерами в рамках данного анализа.  
**Цель:** сохранить OpenMP-выполнение и все режимы NeoQC, добавить гибридный `--compute cuda` для ускорения доказанных вычислительных bottleneck без изменения семантики тримминга и QC.

## 1. Текущая архитектура

### QC-only (single/paired)

1. `FastqReader` читает FASTQ/FASTQ.GZ, формирует batch размером 100 000 записей/пар.
2. `src/main.cpp::processBatchParallel` / `processPairedBatchParallel` распределяет записи по OpenMP threads.
3. У каждого потока свой `QualityAnalyzer`, включая собственную точную карту дупликаций `unordered_map<DuplicationKey,uint64_t>`.
4. После EOF `mergeAnalyzers`: резервирование итоговой карты, сложение метрик, перенос узлов дупликации. Далее сводки, TSV и отчёты.

### С триммингом (`--trim`, `--full-trim`)

1. `processTrimmedBatchParallel()` / `processTrimmedPairedBatchParallel()` выполняют **QC before → trimming → QC after** в одном OpenMP проходе по batch. Не нужно снова перечитывать исходные FASTQ.
2. В PE режиме одновременно поддерживаются **четыре набора** thread-local `QualityAnalyzer`: R1/R2 × before/after.
3. `Trimmer::trimPairDetailed()` применяет UMI, quality trimming, затем (только по явно запрошенным параметрам) overlap/correction, adapter trimming, polyG/polyX, min-length; опционально merge reads. Точный порядок см. в коде `src/trimming/trimmer.cpp`.
4. После OpenMP цикла `FastqWriter` последовательно записывает прошедшие reads в `.fastq.gz`; затем все наборы статистик объединяются и формируется before/after QC.
5. `--full-trim` задаёт профиль адаптеры + quality (Q20, window 4) + polyG/polyX + minimum length 50; **не** включает по умолчанию overlap correction, merge или UMI.
6. `--skip-adapters` относится к adapter-content QC, а **не** к отключению adapter trimming внутри `--full-trim`.

Ключевые файлы:
- `src/main.cpp` — CLI, batch OpenMP, тримминг+QC, merge, вывод.
- `src/quality_analyzer.cpp`, `include/quality_analyzer.h` — QC, точная 50-nt duplication, thread-local stats.
- `src/trimming/trimmer.cpp`, `src/trimming/adapter_trimmer.cpp` — оркестрация тримминга и поиск адаптеров.
- `src/trimming/pe_overlap.cpp`, `src/trimming/pe_processor.cpp` — overlap/correction/merge.
- `src/trimming/quality_trimmer.cpp`, `poly_tail_trimmer.cpp`, `umi_processor.cpp` — другие операции.
- `src/trimming/fastq_writer.cpp` — gzip output.
- `docs/problems/NeoQC_trim_optimization.md` — предыдущая оценка узких мест, требует нового baseline на текущем коммите.

## 2. Правила совместимости

- Два вычислительных режима: `--compute omp` (по умолчанию, старое поведение) и `--compute cuda` (гибрид OpenMP+GPU). Третий serial CPU backend не нужен.
- CUDA необязательна при сборке; OMP-бинарник должен работать без NVIDIA GPU и CUDA Toolkit.
- Все параметры `--trim`, `--full-trim`, quality/adapter/overlap/UMI, `--merge`, `--samples`, `--plot`, `--skip-adapters`, `--timing` работают в обоих режимах.
- `--compute` выбирает **место вычисления**, но не меняет метрики, trimming policy, FASTQ output, порядок reads или формат отчётов.
- Выбор GPU-подзадач не должен допускать двойного подсчёта метрик. Перед первым kernel нужен контракт частичных результатов.
- Точная 50-nt duplication должна остаться точной (без sampling и без скрытой смены ключа).
- В проекте **не добавлять новую метрику k-mer counting**. Уже используемые внутренние k-mer seed/adapter-matching структуры не являются новой QC-метрикой.

## 3. Матрица компонентов: где считать

| Операция | `--compute omp` | Цель для `--compute cuda` | Приоритет и причина |
|---|---|---|---|
| CLI, configuration, samples.csv | CPU | CPU | Оставить |
| FASTQ reading, gzip input, PE-sync | CPU | CPU | Оставить; I/O bottleneck измерять |
| FASTQ validation | CPU/OMP | CPU/OMP на первом этапе | Сохранить детерминированные ошибки |
| Batch creation, staging, output lifecycle | CPU | CPU | Оставить; позже bounded pipeline |
| QC before: per-base quality/histogram | OMP | CUDA | Кандидат 1 после baseline |
| QC after: per-base quality/histogram | OMP | CUDA | То же, только passed reads |
| Per-base A/C/G/T/N, readsPerPosition | OMP | CUDA | Вычислять совместно с quality при возможности |
| Per-sequence quality (rounded/truncated) | OMP | CUDA | После per-base parity |
| Per-sequence GC, GC FastQC-compatible | OMP | CUDA, только при точной parity | Семантика округления/дробных bins критична |
| Read lengths, Q20/Q30, total counts | OMP | CUDA либо OMP — по профилю | Считать вместе с основным kernel, если выгодно |
| Adapter-content QC (before/after) | OMP | Сначала OMP, CUDA только при выгоде | Не путать с удалением адаптеров |
| Adapter trimming (`AdapterTrimmer`) | OMP | Сначала оптимизация OMP; CUDA — отдельный эксперимент | Поиск по позициям и tie-break |
| Quality trimming, sliding window | OMP | Возможно CUDA во второй волне | Требует записи trim boundaries/flags |
| PolyG/PolyX trimming | OMP | Возможно fused CUDA trimming | Не портить последовательность операций |
| UMI extraction | OMP | OMP | Изменения IDs/metadata, слабый кандидат |
| PE overlap detection | OMP | Оставить OMP сначала; возможный CUDA candidate | Только opt-in correction/merge; exact parity |
| Overlap correction и PE merge | OMP | OMP сначала | Структурные изменения пары и качества |
| Min-length filtering, paired retention | OMP | CPU/OMP как контракт результата | После всех trim-шагов |
| Duplication (50 nt exact), overrepresented | OMP | Оставить CPU/OMP до решения memory problem | GPU не решает автоматически OOM; hash-таблицы |
| Per-thread/local results merging | CPU | CPU + компактные GPU reductions | Снизить peak RSS |
| gzip output trimmed FASTQ | CPU | CPU | Оптимизировать буферизацию/конвейер отдельно |
| TSV, QC evaluation, HTML/plots | CPU/Python | CPU/Python | Не переносить |

Рекомендуемая **первая реализация CUDA**: метрики QC (до/после тримминга) на подготовленных read batches. **Не начинать** с GPU-адаптеров или GPU-overlap: нужно сначала доказать, что это лучший узкий участок, и избежать дублирования CPU/GPU-обработки.

## 4. Основной риск — RAM/OOM

Код `src/main.cpp::mergeAnalyzers()` сначала суммирует `getDuplicationEntryCount()`, вызывает `reserveDuplicationEntries(sum)`, затем делает `destination.merge(source)` и `destination.mergeDuplicationCounts(source)`. В `src/quality_analyzer.cpp` использован `unordered_map::extract` / `insert` с переносом узлов. Это **уже не наивное полное копирование** всех элементов, но возможно высокое пиковое выделение bucket-array при `reserve(sum)`, особенно когда между потоками много совпадающих ключей. Кроме того, четыре независимых карты `before/after × R1/R2` могут быть большими.

`Killed` после `Processed ...` **не доказывает** OOM и не доказывает, что виноват `merge`. Проверить OOM killer, peak RSS, фазу и окружение (cgroup/container limit); не рекомендовать 'GPU решит память' без измерений.

Проверка на реальном проблемном образце:

```bash
/usr/bin/time -v -o baseline.time.txt \
  env OMP_NUM_THREADS=20 ./build/neoqc \
  --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> \
  --sample-id baseline --out /path/to/baseline \
  --full-trim --timing --plot

echo "exit=$?"
journalctl -k --since '30 minutes ago' | grep -iE 'oom|killed process|out of memory'
# Если запущено в cgroup v2:
cat /sys/fs/cgroup/memory.events 2>/dev/null
```

Число потоков и доступная RAM для реального теста должны соответствовать проблемной машине. Если падение воспроизводится, отдельно проверить 4, 8, 20 потоков и степень уникальности дупликационных префиксов. Для проверки OOM **важно** иметь реальные kernel/cgroup logs, а не только exit code `137`.

### Инструментация P0

Логировать по `std::chrono::steady_clock` и текущему/пиковому RSS:
- FASTQ read/decode;
- QC before;
- `Trimmer::trimPairDetailed` (внутри — adapter/quality/overlap/poly/UMI по отдельному sampling-трассированию, не на каждой записи);
- QC after;
- writer + gzip;
- `mergeAnalyzers` отдельно для R1-before, R2-before, R1-after, R2-after;
- построение duplication/overrepresented;
- output/report.

Держать замеры малозатратными: время по batch, sampling внутри hot loops, не `clock_gettime()` на каждый read.

## 5. Этапы с конкретным первым результатом

| Этап | Суть | Результат/acceptance gate | Зависимость |
|---|---|---|---|
| **0A — baseline** | Закрепить коммит, 4 режима (QC-only; `--trim`; `--full-trim`; `--full-trim --skip-adapters`), фиксировать cores/RAM/reads | Файл `benchmarks/baseline.md` с wall-time, RSS, thread count, output size, status | Нет |
| **0B — OOM/merge** | Фазовые timings + RSS; диагностировать `Killed`, четыре duplication maps; проверить `reserve(sum)` | Доказана/отвергнута гипотеза OOM, memory profile, регрессионные тесты | 0A |
| **0C — trim bottlenecks** | Замерить adapter matching, overlap(opt-in), gzip output; оптимизация без изменения семантики | Сравнение до/после на trimmed FASTQ и `trimming_report.json` | 0A |
| **1 — CLI/build** | `--compute omp` по умолчанию, optional CMake CUDA; cuda без toolkit даёт понятную ошибку | CPU-only tests проходят, прежняя CLI не ломается | 0A |
| **2 — QC contracts** | Вынести компактные промежуточные QC counters и согласовать before/after, R1/R2, filtered pairs | CPU reference tests для partial stats и merge; без изменения форматов TSV | 1 |
| **3 — CUDA MVP** | Packed batch (offset+length+sequence+quality), H2D, первый kernel пер-base quality+N, compact D2H | CUDA vs OMP exact integer equality на synthetic FASTQ | 2 |
| **4 — CUDA QC before/after** | Два набора метрик без повторного чтения; решить хранение original/trimmed batch или GPU side boundary | Все before/after R1/R2 parity tests; peak RSS ограничен | 3 |
| **5 — остальные QC** | Per-sequence quality, GC (в том числе FastQC distribution), positions, Q20/Q30 | Численная parity с явно заданной точностью float | 4 |
| **6 — trimming CUDA feasibility** | Измерить cost adapter/quality/overlap и H2D/D2H; эксперимент с одним hot kernel, не автоматический port | Решение: оставить OMP либо GPU; end-to-end benchmark | 0C, 4 |
| **7 — pipeline/I/O** | Bounded reader/worker/writer, gzip bulk-buffering и сохранение порядка | Identical decompressed trimmed FASTQ, peak RSS не выше бюджета | 0B, 4 |
| **8 — benchmark/stabilize** | Single/paired, trim/no-trim, large/high-unique, batch modes; fail-fast при отсутствии GPU | OMP и CUDA сборки, regression suite и реальные speedups | 3–7 |

Этапы 0A–0C желательно делать в отдельных небольших коммитах/PR. Не переписывать `main.cpp` целиком в одном PR.

## 6. Первая задача, которую стоит делать прямо сейчас

### Issue 0A: `Add phase timing and RSS instrumentation for trim pipeline`

**Изменить:** `src/main.cpp` (первый шаг), опционально новый `include/perf_metrics.h`/`src/perf_metrics.cpp`.

**Подзадачи:**
1. Расширить `--timing`: read, OpenMP batch, write и post-EOF finalize без изменения обычного вывода.
2. Внутри OpenMP собирать времена QC-before, trim, QC-after как sampled counters; отдельно общий wall time batch, чтобы не складывать thread CPU-times как wall time.
3. При входе/выходе четырёх `mergeAnalyzers` логировать RSS и число duplication keys.
4. Поддержать CSV/TSV benchmark log или компактный текстовый summary по окончании, даже в OMP режиме.
5. Тест: обычные QC результаты и outputs должны быть идентичны запуску без timing; `ctest --test-dir build --output-on-failure`.
6. Прогнать 1 небольшой fixture и 1 реальный проблемный dataset; записать аппаратную конфигурацию, команду, `time -v`, OOM/cgroup logs.

**Definition of done:** известна фаза максимума RSS и процент wall time, приходящийся на `trim`, `QC`, gzip write и finalize. После этого выбирается узкий участок. Никакой CUDA в этой задаче нет.

### Issue 0B: `Reduce exact duplication merge peak memory without changing semantics`

Проводить **после** измерений. Проверить эффект `reserve(sum)` на высокой доле межпоточных повторов, переполнение/избыточную ёмкость, память четырёх наборов, а также стоимость `rehash(0)`. Рассмотреть перенос карт по одной и резервацию с учётом фактических distinct (без отключения exact metric). Тестировать полное совпадение duplication TSV/JSON.

### Issue 1: `Introduce optional OMP/CUDA backend selection`

После стабильного baseline и memory диагностики. Добавить `--compute omp|cuda`, optional CUDA compile; при режиме `omp` старый поток вычислений нетронут. Не делать фиктивный CUDA backend, который всё равно считает всё на CPU без сообщения пользователю.

## 7. CUDA архитектура, минимально необходимая

```
Reader (CPU) -> paired ReadBatch (bounded)
  |                        |
  |    --compute omp       |    --compute cuda (hybrid)
  v                        v
OpenMP: before QC       CUDA: before QC compact counters
OpenMP: trim            OpenMP: trim (at first)
OpenMP: after QC        CUDA: after QC on passing pairs
  |                        |
  +-------- partial QC stats + exact duplication CPU maps ------+
                               |
                      compact merge and reports
```

Не дублировать GPU QC и `QualityAnalyzer::processRecord()` для одних и тех же счётчиков. Для duplication можно предусмотреть отдельный CPU-проход по batch/коротким ключам, а GPU использовать только для выбранных QC counters. Это требует декомпозиции `processRecord` на независимые части. Общий интерфейс должен работать на batch и возвращать **частичные статистики**; не требуется немедленно оборачивать весь `QualityAnalyzer` в virtual class.

Для GPU данные следует упаковать в SoA/flat buffers: sequence bytes, qualities, offsets, lengths, flags passed, mate IDs, исходные/итоговые границы. Нельзя пересылать `std::string`/`unordered_map` на GPU как есть. Сначала предпочесть один host copy и один или несколько kernels; пиковую RAM и VRAM контролировать, не удерживая две лишние полные копии больших batch.

Варианты архитектуры before/after:
- **MVP:** копирование двух packed representations для QC-before и QC-after после CPU-тримминга; просто, но дополнительная копия и больший peak RSS. Проверить на ограниченном batch.
- **Следующий этап:** GPU-тримминг, возвращающий границы/flags, только если end-to-end benchmark подтверждает эффективность и полное соблюдение контракта.

## 8. Приоритеты и запреты

1. **Сейчас:** замеры/OOM/merge → устранение измеренного CPU bottleneck тримминга или memory problem.
2. **Затем:** optional backend + CUDA QC MVP и before/after с доказанной parity.
3. **Только по результатам профилирования:** оптимизация gzip pipeline, GPU adapter trimming, GPU overlap.
4. **Не делать:** K-mer counting metric, GPU I/O ради GPU, скрытую подмену точной duplication approximation, отключение before/after QC ради скорости, сравнение только CUDA kernel-time вместо полного wall time.

## 9. Проверочная матрица

Для всех существенных изменений обязательны single-end/paired-end, clean/dirty, mixed-length/short, адаптерные positive/negative, N/GC/Phred boundaries, `--trim`, `--full-trim`, `--overlap-correction`, `--merge`, UMI, `--skip-adapters`, `--plot`, `--samples`. Для QC — integer counts **exact match**, FP — явная tolerance по согласованной семантике, FASTQ — полное совпадение sequence+quality+order (при сравнении gzip сначала распаковать), reports — неизменные значения и структура.

Исходная ветка: https://github.com/kk2gvb/NeoQC/tree/optimization_and_cuda
