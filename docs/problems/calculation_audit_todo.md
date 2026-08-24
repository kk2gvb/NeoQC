# NeoQC: аудит расчётной части и TODO

Дата аудита: 2026-08-24.

Область проверки: чтение FASTQ, базовые счётчики, Phred-метрики,
GC/base/N content, длины ридов, адаптеры, точная дупликация, OpenMP merge,
формирование TSV и вычисление статусов по ruleset
`fastqc-compatible-v1`.

Этот документ фиксирует проблемы и необходимые проверки. Расчётный код в
рамках данного аудита не менялся. Риск памяти exact duplication и будущая
RNA-seq policy подробно вынесены в соседний документ
[`duplication_memory_and_rna_qc_policy_todo.md`](duplication_memory_and_rna_qc_policy_todo.md).

## P0 — целостность результата

### 1. Повторный запуск смешивает новые и старые артефакты

Подтверждено воспроизведением: после paired-end запуска и последующего
single-end запуска в том же `--out` остаются все TSV и summary для R2.
Аналогично, запуск с `--skip-adapters` не удаляет ранее созданные
`adapter_content_R1.tsv` и `adapter_content_R2.tsv`.

Это влияет на расчёт статуса, а не только на внешний вид каталога:
`scripts/qc_rules.py::_active_reads()` считает направление активным, если
находит хотя бы один известный TSV. Поэтому новый R1 может быть оценён вместе
со старым R2. Старый adapter TSV также будет оценён отдельным запуском
`evaluate_qc.py`, хотя в новом запуске адаптеры не вычислялись.

Риск распространяется и на режим `--plot`: запуск без него не удаляет старые
`qc_evaluation.json`, `plots_manifest.json`, графики и компактный NeoQC report.
Большинство TSV пишется напрямую через `ofstream`; при ошибке записи возможен
каталог из файлов разных поколений. Транзакционный marker и atomic publication
сейчас реализованы только для duplication artifacts.

TODO:

1. Ввести manifest одного запуска с `run_id`, направлением (`R1`/`R2`),
   параметрами и перечнем успешно опубликованных файлов.
2. Перед публикацией результата удалить или изолировать артефакты направлений
   и метрик, отсутствующих в новом запуске.
3. Писать весь набор результатов во временный каталог и атомарно публиковать
   его только после успешного окончания расчёта.
4. Научить rule-engine определять активные reads по manifest, а не по любому
   случайно найденному TSV.
5. Добавить regression-сценарии paired -> single, adapters -> skip-adapters,
   plot -> no-plot и успешный запуск -> аварийный повторный запуск.

Критерий готовности: один каталог никогда не может содержать логически
активную смесь двух запусков; оценка нового single-end результата не содержит
R2.

### 2. Batch `case_summary.json` не обновляется

В batch-режиме summary создаётся только при отсутствии файла:
`if (!fs::exists(caseOutDir / "case_summary.json"))`. Повторный расчёт того же
пациента оставляет старые статусы, числа ридов, ошибки и дату. Кроме того,
`neoqc_version` в этом файле жёстко задан как `0.1`, тогда как актуальная версия
проекта — `1.0.1`.

TODO:

1. Всегда формировать новый case summary из результатов текущего batch run.
2. Публиковать его атомарно после окончания всех sample runs.
3. Получать версию из единственного build-time источника, а не из строкового
   литерала.
4. Записать в summary ruleset ID, ruleset version/SHA-256, параметры запуска и
   `run_id`.
5. Добавить тест повторного batch-запуска с изменившимся результатом sample.

### 3. FASTQ quality encoding не валидируется

Расчёт безусловно использует `quality_byte - 33`, но reader проверяет только
непустую строку и совпадение её длины с sequence. Нет проверки допустимого
ASCII-диапазона и нет явного отказа для несовместимой кодировки.

Следствия:

- Phred+64 будет молча интерпретирован как завышенный Phred+33;
- байт с вычисленным `q > 93` попадает в `qualitySum`, `qualityCount`, Q20/Q30
  и среднее по риду, но не попадает в histogram размером 94;
- quantile затем использует полный `qualityCount`, не достигает target по
  histogram и возвращает последний bin 93;
- байт с `q < 0` молча исключается из quality statistics, тогда как знаменатель
  Q20/Q30 остаётся равным числу оснований.

В итоге разные quality-метрики одного и того же файла могут быть внутренне
несогласованными.

TODO:

1. Зафиксировать поддерживаемую кодировку: Sanger/Illumina Phred+33.
2. До расчёта отклонять каждый quality byte вне поддерживаемого printable
   диапазона с номером записи и позиции.
3. Если нужна поддержка старого Phred+64 — сделать её только через явный
   параметр либо надёжно протестированное определение кодировки, сохраняя
   выбранную кодировку в provenance.
4. Добавить тесты на `!` (Q0), `~` (Q93), управляющий/DEL byte и Phred+64
   fixture.
5. Проверять инварианты: сумма histogram равна `qualityCount`, сумма
   per-sequence distribution равна total reads.

## P1 — расхождения расчётов и профиля `fastqc-compatible-v1`

### 4. N входит в знаменатель base sequence content

`per_base_sequence_content_R*.tsv` делит A, C, G и T на
`A + C + G + T + N`. Официальный модуль FastQC нормирует A/C/G/T только на
сумму канонических оснований, а N рассматривает отдельной метрикой.

Из-за этого N-rich позиция искусственно уменьшает `abs(A-T)` и `abs(G-C)`.
Например, при 10% A и 90% N NeoQC получает разницу A/T 10 процентных пунктов,
тогда как среди известных оснований она равна 100. Это может скрыть WARN/FAIL
по `maximum_base_difference_percent`.

Официальная реализация для сопоставления:
[FastQC PerBaseSequenceContent.java](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/PerBaseSequenceContent.java).

TODO:

1. Развести два знаменателя: A/C/G/T — среди canonical calls, N — среди всех
   ридов, покрывающих позицию.
2. Определить поведение позиции, содержащей только N: `NOT_EVALUATED`/NaN либо
   явно документированный neutral value.
3. Добавить golden tests с 0%, 50%, 90% и 100% N и сравнение с выбранной
   версией FastQC.

### 5. Наблюдаемое GC-распределение не эквивалентно FastQC

NeoQC присваивает каждому риду ровно один bin через
`lround(GC * 100 / length)`. FastQC использует модель для конкретной длины
рида и дробно распределяет наблюдение по соседним процентным bins. Python-слой
NeoQC строит близкую теоретическую normal curve уже поверх другого observed
distribution.

Расхождение особенно заметно для коротких и переменной длины ридов и может
изменить `modeled_deviation_percent`, следовательно WARN/FAIL.

Официальная реализация для сопоставления:
[FastQC PerSequenceGCContent.java](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/PerSequenceGCContent.java).

TODO:

1. Зафиксировать конкретную референсную версию FastQC для compatibility
   profile.
2. Перенести либо независимо воспроизвести length-aware fractional GC model.
3. Создать golden corpus для коротких, длинных и mixed-length reads и сравнить
   observed bins, theoretical curve, deviation и итоговый статус.
4. Если текущий алгоритм сохраняется намеренно, не называть этот модуль строго
   FastQC-compatible и версионировать новую семантику.

### 6. Per-base quality получает статус на слишком малой выборке

NeoQC вычисляет quartile и median при любом `qualityCount > 0`, отдельно для
каждой позиции. FastQC не публикует percentile для позиции/группы, если в ней
недостаточно наблюдений (в текущем официальном коде требуется более 100), а для
длинных ридов использует группировку позиций.

Поэтому один или несколько ридов могут дать NeoQC WARN/FAIL там, где
FastQC-compatible статистика не должна считаться надёжно оценённой. Отдельные
позиции вместо FastQC grouping также меняют minimum quartile/median.

Официальная реализация для сопоставления:
[FastQC PerBaseQualityScores.java](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/PerBaseQualityScores.java).

TODO:

1. Определить minimum observation count и состояние `NOT_EVALUATED`.
2. Формально решить, нужна ли FastQC BaseGroup semantics или собственная
   per-position semantics.
3. Добавить тесты на 1, 100, 101 read и на mixed-length хвосты с малым
   покрытием.
4. Проверить quartile definition на граничных малых наборах против выбранной
   версии FastQC.

### 7. Набор адаптеров не соответствует официальному default-набору FastQC

NeoQC ищет точный 12-nt prefix каждого встроенного адаптера. При этом:

- нет PolyA и PolyG, присутствующих в актуальном default adapter list FastQC;
- NeoQC SmallRNA5 начинает поиск с `GTTCAGAGTTCT`, тогда как FastQC использует
  `GATCGTCGGACT`;
- версия/контрольная сумма референсного списка нигде не закреплена;
- пользователь не может передать собственный adapter list.

Следовательно, одинаковый FASTQ может получить другое maximum adapter content
и другой статус при профиле с названием `fastqc-compatible-v1`.

Официальные источники для сопоставления:
[adapter_list.txt](https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt),
[AdapterContent.java](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/AdapterContent.java).

TODO:

1. Утвердить встроенный adapter set и его ориентацию отдельно для R1/R2.
2. Добавить PolyA/PolyG, если целью остаётся совместимость с выбранной версией
   FastQC.
3. Разрешить versioned пользовательский adapter config.
4. Добавить положительный и отрицательный fixture на каждый adapter, включая
   несколько вхождений и короткие reads.

### 8. Exact duplication и FastQC duplication — разные метрики

NeoQC считает все уникальные 50-nt префиксы exact. FastQC использует
историческую ограниченную схему наблюдения/экстраполяции. Возвращать лимит
100 000 нельзя: это ухудшит научную семантику NeoQC и противоречит принятому
exact-контракту. Но одинаковые пороги под именем `fastqc-compatible-v1` не
гарантируют одинаковый статус для этих двух разных оценок.

TODO:

1. Сохранить exact как основной алгоритм и явно записывать algorithm ID в
   evaluation/report, не только в duplication summary.
2. Описать профиль как «FastQC-like thresholds over NeoQC exact metric» либо
   выделить отдельный NeoQC ruleset.
3. Подготовить сравнительные fixtures: low complexity, high complexity,
   дубликаты, распределённые дальше первых 100 000 уникальных reads.
4. Не выдавать возможный будущий approximate mode за exact.

Официальная реализация для сопоставления:
[FastQC DuplicationLevel.java](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/DuplicationLevel.java).

## P2 — контракт, границы и полнота валидации

### 9. Overrepresented sequences: несовпадающая граница и отсутствующий статус

Код включает последовательность только при доле строго больше `0.1%`.
`docs/sequence-duplication.md` говорит «above 0.1%», но
`docs/QC_METRICS.md` — «не менее 0.1%» и задаёт WARN начиная с `0.1%`.
Ровно 1 read из 1000 поэтому трактуется по-разному кодом и документацией.

Кроме того, ruleset вообще не содержит правила для overrepresented sequences:
TSV создаётся, но описанные в `docs/QC_METRICS.md` PASS/WARN/FAIL не
вычисляются. Поле `possible_source` всегда равно `No Hit`; сопоставления с
contaminant database нет.

TODO:

1. Утвердить включительность границы и покрыть значения ниже, ровно и выше
   0.1% тестами.
2. Либо добавить версионированное правило статуса, либо удалить обещание
   статусов из документации.
3. Переименовать `possible_source` в нейтральное поле/удалить его до появления
   реального matching либо реализовать versioned contaminant database.

### 10. Length rule содержит недостижимую FAIL-ветку

Reader отклоняет пустую sequence, поэтому native NeoQC не может сформировать
валидный `minimum_length <= 0`. Соответствующий FAIL check ruleset фактически
проверяет только вручную созданный или повреждённый TSV. Одновременно любое
число наблюдаемых длин больше одного даёт WARNING, хотя variable length может
быть ожидаемым результатом trimming или протокола.

TODO:

1. Отделить проверку целостности TSV от научного QC rule.
2. Решить, является ли variable length универсальным warning или
   library/protocol-specific observation.
3. Добавить fixtures с fixed length, несколькими длинами и редким коротким
   хвостом.

### 11. Rule-engine недостаточно проверяет внутреннюю согласованность TSV

Общий parser проверяет числовой формат, конечность и неотрицательность, но не
проверяет предметные инварианты:

- позиции/bins должны быть целыми, уникальными и упорядоченными;
- A+C+G+T+N должны соответствовать определённому контрактом знаменателю;
- counts должны быть целыми, а их суммы — согласованы с total reads;
- cycle coverage не должна расти после более короткой позиции;
- одинаковые метрики R1/R2 должны принадлежать одному run.

Повреждённый или смешанный TSV поэтому может получить формально валидный
PASS/WARN/FAIL.

TODO:

1. Добавить metric-specific schema validation и cross-file invariants.
2. Проверять run identity и source checksum через manifest.
3. Расширить malformed tests на duplicate/missing bins, дробные counts,
   неверные суммы процентов и несовпадающие totals.

### 12. Нет воспроизводимого численного regression-suite из чистого checkout

Все файлы `tests/data/*` игнорируются правилом `.gitignore: data`. CMake
регистрирует много тестов с этими путями, но не запускает
`scripts/generate_test_data.sh` как prerequisite. Локально fixtures могут
существовать, а в чистом checkout соответствующие CTest-тесты падают из-за
отсутствующих входов и не проверяют расчёты.

Также отсутствует единый differential test, сравнивающий все численные TSV при
`OMP_NUM_THREADS=1` и `OMP_NUM_THREADS=N`. В ходе этого аудита существующий
medium paired fixture дал побайтово одинаковые outputs при 1 и 4 потоках, то
есть конкретное расхождение не обнаружено; проблема состоит в отсутствии
постоянной защиты этого инварианта.

TODO:

1. Либо хранить минимальные deterministic fixtures в Git через точное
   исключение из `.gitignore`, либо генерировать их отдельной CMake-командой с
   declared outputs/dependencies.
2. Не использовать генератор, который без отдельной защиты удаляет всё
   содержимое широкого `tests/data`.
3. Добавить небольшие golden fixtures с вручную проверенными ожидаемыми
   значениями для каждой метрики.
4. Добавить differential OMP test для single/paired, разных длин, адаптеров и
   duplication merge.
5. Запускать C++ tests и Python rules/report tests в CI; plot tests включать в
   окружении с зафиксированной версией matplotlib.

## Рекомендуемый порядок исправления

1. Изоляция/атомарная публикация запуска и обновляемый `case_summary.json`.
2. Строгая проверка Phred+33 и внутренние quality invariants.
3. Исправление знаменателя per-base sequence content.
4. Формальное решение о степени FastQC compatibility для GC, per-base quality,
   adapters и exact duplication.
5. Metric-specific validation входных TSV.
6. Воспроизводимый golden/differential test suite из чистого checkout.

До выполнения пунктов 1–3 результат из повторно используемого каталога и
результат на FASTQ с неизвестной quality encoding нельзя считать надёжно
самодостаточным без ручной проверки provenance.
