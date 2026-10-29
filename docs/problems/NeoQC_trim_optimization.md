# NeoQC: оптимизация полного тримминга

Дата: 2026-10-27.

## Цель

Сократить время paired-end запуска с `--full-trim`, сохранив неизменными:

- FASTQ до и после тримминга;
- полный QC-отчёт для исходных и итоговых reads;
- `trimming_report.json`;
- точную проверку пар, правила фильтрации и воспроизводимость результата.

Оптимизация не должна менять публичную семантику CLI и результаты QC без
отдельного согласованного изменения версии формата/алгоритма.

## Наблюдаемая проблема

На образце `CHE_S6` были получены следующие повторяемые результаты (3 прогона,
20 OpenMP threads):

| Режим | Среднее время |
|---|---:|
| NeoQC QC без тримминга | 91.926 s |
| NeoQC `--full-trim` | 1049.467 s (17.49 min) |
| fastp preprocessing | 36.407 s |
| cutadapt preprocessing | 43.348 s |

Вход:

- R1: 44 673 726 reads, 101 bp, 2.5 GiB gzip;
- R2: 44 673 726 reads, 87 bp, 2.3 GiB gzip.

Источники измерений: `benchmark_results/raw_times.csv` и
`benchmark_cutadapt/raw_times.csv`.

Важно: разница между 91.9 s и 1049.5 s не является чистым временем алгоритма
тримминга. В `--full-trim` дополнительно есть QC после тримминга, запись двух
gzip FASTQ и финальная агрегация QC-структур.

## Что не является ошибкой

Полные отчёты *до* и *после* нужны продукту. Их нельзя убрать ради скорости и
некорректно называть самостоятельным устранимым bottleneck.

Для точного отчёта после тримминга данные surviving reads должны быть учтены
отдельно: это минимально означает обработку исходного read и его итогового
среза. Цель оптимизации — не отказаться от этой работы, а:

1. не выполнять дополнительные проходы, аллокации и синхронизации вокруг неё;
2. не хранить и не сортировать больше данных, чем требует точная метрика;
3. перекрыть CPU-обработку с вводом-выводом;
4. резко сократить стоимость поиска overlap и адаптеров.

## Текущая схема

```text
R1.gz + R2.gz
       │
       ├─ последовательное чтение, FASTQ/Phred/pair validation
       │
       ├─ batch 100 000 пар
       │    │
       │    ├─ OpenMP: QC before R1/R2
       │    ├─ OpenMP: trim pair
       │    └─ OpenMP: QC after R1/R2 для passed pairs
       │
       └─ последовательная запись R1.gz + R2.gz

после EOF: merge exact-duplication maps → sort → merge → отчёты
```

QC before/after в текущем виде встроен в один OpenMP проход по batch, что
правильно: отдельного повторного чтения исходного FASTQ для after-QC нет.
Однако `QualityAnalyzer` всё равно проходит по bases исходного и итогового
read отдельно, а exact duplication собирается для четырёх потоков данных:
R1/R2 × before/after.

## Узкие места

### P0. Paired-end overlap search для каждой пары

До исправления adapter trimming неявно вызывал `PeOverlapAnalyzer::find()` для
каждой пары. Алгоритм перебирает каждый допустимый shift и посимвольно
сравнивает overlap.

Для измеренного образца (R1=101, R2=87, `min_overlap=30`) это 129 кандидатов
на пару и до 7 917 сравнений bases на пару в худшем случае. Дополнительно на
каждую пару создаётся reverse complement R2.

Риск: десятки миллионов пар превращают даже короткие посимвольные циклы в
основную CPU-нагрузку.

Решение: обычный adapter trimming теперь не запускает поиск перекрытия.
`PeOverlapAnalyzer::find()` вызывается только при явно заданных
`--overlap-correction` или `--merge`; read-through trimming через overlap не
входит в `--full-trim`.

Код:

- `src/trimming/trimmer.cpp:81-85`;
- `src/trimming/pe_overlap.cpp:80-132`.

### P0. Adapter search имеет полный перебор и аллокации на read

Если read-through overlap не найден, `AdapterTrimmer` для каждого R1 и R2:

1. создаёт `std::vector<std::string>`;
2. копирует в него все adapter sequences;
3. проверяет каждую позицию каждого адаптера;
4. не прекращает поиск после первого найденного кандидата, поскольку выбирает
   глобально лучший match.

В standard profile семь адаптеров. Этот путь выполняется для огромной доли
reads, где адаптера нет.

Текущий код: `src/trimming/adapter_trimmer.cpp:53-99`.

### P1. Сериализованные gzip I/O и запись по одному read

Чтение R1/R2, проверка идентификаторов и запись выполняются одним потоком.
После каждого OpenMP batch основной поток последовательно сериализует passed
records в строку и вызывает `gzwrite` для каждого read.

Это создаёт барьер вида «сначала прочитать → затем вычислить → затем записать»:
ввод, compute и compression не перекрываются. В output path также происходит
несколько динамических конкатенаций строк на read.

Текущий код:

- `src/main.cpp:1615-1638`;
- `src/trimming/fastq_writer.cpp:25-55`.

### P1. Exact duplication: хеш-таблицы, сортировка и merge

Для отчётов before/after требуется точная duplication-метрика. Сейчас для
каждого OpenMP worker поддерживается отдельная `unordered_map`, затем записи
копируются, сортируются и сливаются. В full-trim это происходит для R1/R2 и
before/after.

Это функционально обоснованная работа, но её memory allocation и финальная
сортировка могут стать заметным хвостом запуска, особенно при высокой
уникальности reads.

Текущий код:

- `src/quality_analyzer.cpp:256-266`;
- `src/main.cpp:1398-1405`.

### P2. Последовательная FASTQ validation и лишние сканирования строк

До OpenMP reader проверяет каждый quality byte. Затем QC снова обрабатывает
quality при расчёте метрик. `readLine()` также копирует буфер в `std::string`
и повторно ищет newline во всей строке. Это не объясняет 17 минут само по
себе, но масштабируется с числом bases и ограничивает дальнейшее ускорение.

Текущий код: `src/fastq_reader.cpp:26-68`, `src/fastq_reader.cpp:180-202`.

### P2. Adapter-content QC

Отдельно от удаления адаптеров QC ищет adapter k-mer в original и trimmed
reads. Это дополнительная работа, но не основной источник разрыва.

`--skip-adapters` отключает только эту QC-метрику. Он **не** выключает
удаление адаптеров в `--full-trim` и потому не устраняет P0 алгоритмы выше.

## Предлагаемая стратегия

### Этап 1. Ускорить adapter matching

1. Хранить готовый список `string_view`/ссылок на адаптеры в
   `AdapterTrimmer`, не копировать его на read.
2. Сначала проверять редкие и длинные seed/k-mer адаптера, затем подтверждать
   только кандидата; не запускать full comparison на каждой позиции.
3. Использовать precomputed uppercase/2-bit representation либо таблицу
   нормализации bases вместо `std::toupper` в горячем цикле.
4. Сохранить текущее правило выбора лучшего match и добавить golden/regression
   tests на совпадающие позиции, mismatch limits и tie-break.

Критерий готовности: идентичные trimmed FASTQ и `trimming_report.json` на
fixture-наборе; отдельный microbenchmark adapter-negative reads быстрее
текущего baseline.

### Этап 2. Ускорить PE overlap

Нужно сначала зафиксировать требуемую семантику: применяется ли overlap только
для read-through adapter trim, либо также для correction/merge.

Этот этап относится только к явно запрошенным correction/merge. Для этих
режимов можно отсеивать большинство пар дешёвым seed-фильтром и подтверждать
только немногочисленные кандидаты. Полный exhaustive поиск следует оставить
как reference implementation для тестов.

Критерий готовности: совпадение найденного `PeOverlap` с reference на
положительных, отрицательных и пограничных fixtures; измеримое сокращение
общего времени `--full-trim` в benchmark matrix.

### Этап 3. Сохранить before/after QC, но снизить его накладные расходы

Отчёты не отключаются. Возможные безопасные улучшения:

1. Сохранить текущий single-pass-per-batch порядок: before QC → trim → after
   QC. Не вводить второй проход по файлу.
2. Вынести structural validation и quality-range validation в общий
   посимвольный цикл QC там, где это не ухудшает раннюю диагностику. Ошибки
   должны сохранять номер файла, read и позиции.
3. Уменьшить число независимых записей в per-base счётчики и использовать
   cache-friendly layout; рассмотреть локальные компактные counters с
   reduction после batch.
4. Исследовать вычисление after-метрик через вычитание удалённых segment-ов
   из before-агрегатов только для метрик, для которых это строго эквивалентно.
   Для duplication, per-sequence GC и quality distributions это не замена
   after-QC и требует отдельного доказательства корректности.

Таким образом, пункт не состоит в «убрать QC after»: он ограничивает лишнюю
работу вокруг обязательных двух наборов статистик.

### Этап 4. Развязать pipeline и gzip

Сделать bounded producer/consumer pipeline с сохранением порядка reads:

```text
reader/decompress → bounded batches → OpenMP trim+QC → ordered writer/compress
```

Writer получает уже сериализованные крупные buffers, а не один `gzwrite` на
record. Число in-flight batches ограничивается, чтобы не ухудшить peak RSS.
Порядок R1/R2 и output должен оставаться детерминированным.

Критерий готовности: byte-level FASTQ equivalence после decompression,
стабильный порядок records, ограниченная память и сокращение общего времени
`--full-trim`.

### Этап 5. Отдельно оптимизировать exact duplication

Не заменять exact duplication sampling-ом без продуктового решения. Сначала
измерить число unique prefixes, RSS и время финальной сортировки. Возможные
направления: reserve по batch, memory resource/arena, по-thread sorted runs и
k-way merge без промежуточного полного копирования.

## Порядок реализации

1. Зафиксировать benchmark matrix и текущий общий baseline.
2. Этап 1: adapter matcher — наименьший риск и наиболее явные лишние
   аллокации.
3. Этап 2: overlap search — ожидаемый основной CPU выигрыш.
4. Этап 4: buffered I/O pipeline.
5. Этапы 3 и 5 — после подтверждения общим benchmark, с обязательными
   regression и memory-тестами.

## Benchmark matrix

Каждую реализацию проверять минимум в режимах:

| Режим | Задача |
|---|---|
| QC only | базовая стоимость QC и duplication |
| `--trim` без trim-операций | QC before/after + gzip output |
| `--full-trim --skip-adapters` | алгоритмы trim без Adapter-content QC |
| `--full-trim` | реальный end-to-end сценарий |

Для каждого режима: 3 прогона, фиксированные `OMP_NUM_THREADS`, cold/warm
cache отмечать отдельно и фиксировать общий wall time, CPU time, peak RSS и
input/output bytes. Функциональный gate: тесты проекта и сравнение
decompressed output/JSON reports с reference для набора fixtures.

## Не делать без отдельного решения

- не отключать before/after QC;
- не менять exact duplication на sampling;
- не менять adapter/overlap semantics ради скорости без versioned contract;
- не включать `--skip-adapters` как скрытую часть `--full-trim`;
- не объявлять ускорение по одному warm-cache прогону.
