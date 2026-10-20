"""Russian translations for the NeoQC HTML report and chart labels.

English is the source language of every string produced by NeoQC. The report
renders both languages; the viewer switches between them without regeneration.
Unknown strings (for example, messages from a custom ruleset) are shown in
English in both modes.
"""
from __future__ import annotations

import re

UI = {
    "Report sections": "Разделы отчёта",
    "Order": "Порядок",
    "Triage": "Приоритет",
    "Sequential": "По порядку",
    "Collapse all": "Свернуть всё",
    "Expand all": "Развернуть всё",
    "Print / Save PDF": "Печать / PDF",
    "Sample": "Образец",
    "Reads": "Риды",
    "Generated": "Создан",
    "Ruleset": "Набор правил",
    "Sequencing quality control": "Контроль качества секвенирования",
    "Read sets": "Наборы ридов",
    "Charts available": "Доступно графиков",
    "Rendering errors": "Ошибки построения",
    "None": "Нет",
    "Technical QC overview · PASS / WARNING / FAIL distribution": "Технический QC · распределение PASS / WARNING / FAIL",
    "QC overview": "Обзор QC",
    "QC module": "Модуль QC",
    "Basic Statistics": "Основная статистика",
    "Basic statistics: before / after": "Основная статистика: до / после",
    "Trimming overview": "Обзор тримминга",
    "Before / after trimming": "До / после тримминга",
    "Before trimming": "До тримминга",
    "After trimming": "После тримминга",
    "Before": "До",
    "After": "После",
    "No result": "Нет результата",
    "Change": "Изменение",
    "Status": "Статус",
    "All reads": "Все риды",
    "Trimming report": "Отчёт о тримминге",
    "Trimming parameters": "Параметры тримминга",
    "Reads and bases": "Риды и основания",
    "Adapters": "Адаптеры",
    "Input reads": "Ридов на входе",
    "Kept reads": "Сохранено ридов",
    "Discarded reads": "Отброшено ридов",
    "Bases removed": "Удалено оснований",
    "Mean length R1": "Средняя длина R1",
    "Mean length R2": "Средняя длина R2",
    "QC status before and after trimming": "Статусы QC до и после тримминга",
    "Trimming removes adapters, low-quality ends and too short reads. QC after trimming is calculated on the trimmed FASTQ with the same ruleset; a status change shows the effect of trimming.": "Тримминг удаляет адаптеры, концы низкого качества и слишком короткие риды. QC после тримминга рассчитывается по обрезанным FASTQ с тем же набором правил; изменение статуса показывает эффект тримминга.",
    "Trimming parameters were not recorded for this run.": "Параметры тримминга для этого запуска не записаны.",
    "Library layout": "Тип библиотеки",
    "Paired": "Парная",
    "Single": "Одиночная",
    "Yes": "Да",
    "No": "Нет",
    "Fixed trim at 5′ end (bp)": "Фиксированная обрезка 5′-конца (п.н.)",
    "Fixed trim at 3′ end (bp)": "Фиксированная обрезка 3′-конца (п.н.)",
    "Quality trim from 5′ end": "Обрезка по качеству с 5′-конца",
    "Quality trim from 3′ end": "Обрезка по качеству с 3′-конца",
    "Sliding-window quality trim": "Обрезка по качеству скользящим окном",
    "Quality threshold (Phred)": "Порог качества (Phred)",
    "Window size (bp)": "Размер окна (п.н.)",
    "Adapter trimming": "Удаление адаптеров",
    "Adapter sequences configured": "Настроено последовательностей адаптеров",
    "Minimum adapter match (bp)": "Минимальное совпадение с адаптером (п.н.)",
    "Minimum read length (bp)": "Минимальная длина рида (п.н.)",
    "Terminal polyG trimming": "Удаление terminal polyG",
    "Terminal polyX trimming": "Удаление terminal polyX",
    "Minimum polyG length (bp)": "Минимальная длина polyG (п.н.)",
    "Minimum polyX length (bp)": "Минимальная длина polyX (п.н.)",
    "PE overlap correction": "Коррекция перекрытия парных ридов",
    "Merge overlapping pairs": "Объединение перекрывающихся пар",
    "5′ UMI extraction": "Извлечение 5′ UMI",
    "UMI length (bp)": "Длина UMI (п.н.)",
    "Trimmed R1 FASTQ.gz": "Обрезанный R1 FASTQ.gz",
    "Trimmed R2 FASTQ.gz": "Обрезанный R2 FASTQ.gz",
    "Merged FASTQ.gz": "Объединённый FASTQ.gz",
    "Trimming statistics were not found.": "Статистика тримминга не найдена.",
    "Reads before": "Ридов до",
    "Reads after": "Ридов после",
    "Too short after trimming": "Короче порога после тримминга",
    "Adapter-trimmed reads": "Риды с удалённым адаптером",
    "Quality-trimmed reads": "Риды с обрезкой по качеству",
    "PolyG-trimmed reads": "Риды с удалённым polyG",
    "PolyX-trimmed reads": "Риды с удалённым polyX",
    "Bases before": "Оснований до",
    "Bases after": "Оснований после",
    "Corrected reads": "Скорректировано ридов",
    "Corrected bases": "Скорректировано оснований",
    "Merged pairs": "Объединено пар",
    "Output reads": "Ридов на выходе",
    "Discarded pairs": "Отброшено пар",
    "No adapter positions were recorded.": "Позиции адаптеров не записаны.",
    "Adapter start position (0-based)": "Позиция начала адаптера (с 0)",
    "Trimmed reads": "Обрезанные риды",
    "No comparison charts are available.": "Сравнительные графики недоступны.",
    "Metric": "Метрика",
    "Needs attention": "Требует внимания",
    "Passed": "Пройдено",
    "Passed · folded": "Пройдено · свёрнуто",
    "Sequence screen · Overrepresented sequences": "Скрининг · Сверхпредставленные последовательности",
    "Sequence": "Последовательность",
    "Count": "Количество",
    "Percentage": "Доля",
    "No sequences exceeded the reporting threshold.": "Ни одна последовательность не превысила порог.",
    "Generated by NeoQC": "Создано NeoQC",
    "Theme": "Тема",
    "Light": "Светлая",
    "Dark": "Тёмная",
    "Language": "Язык",
    "Close": "Закрыть",
    "Summary files were not found.": "Файлы сводки не найдены.",
    "Adapter analysis was disabled for this run.": "Поиск адаптеров был отключён в этом запуске.",
    "The source TSV was not produced.": "Исходный TSV не был создан.",
    "Chart was not generated.": "График не построен.",
    "Chart asset unavailable": "Файл графика недоступен",
    "QC evaluation is not available.": "Оценка QC недоступна.",
    "PASS / WARNING / FAIL are technical QC flags from the displayed versioned ruleset, not clinical conclusions. "
    "Plot availability errors are reported separately and never converted into biological FAIL results.":
        "PASS / WARNING / FAIL — технические флаги QC по указанному версионированному набору правил, а не клиническое "
        "заключение. Ошибки построения графиков указываются отдельно и никогда не превращаются в биологический FAIL.",
    # module titles
    "Per base sequence quality": "Качество по позициям",
    "Adapter content": "Содержание адаптеров",
    "Per base sequence content": "Нуклеотидный состав по позициям",
    "Per sequence GC content": "GC-состав ридов",
    "Per base N content": "Доля N по позициям",
    "Sequence length distribution": "Распределение длин ридов",
    "Sequence duplication levels": "Уровни дупликации",
    "Per sequence quality scores": "Среднее качество ридов",
    # basic statistics
    "Processed reads": "Обработано ридов",
    "Total bases": "Всего оснований",
    "Min length": "Мин. длина",
    "Max length": "Макс. длина",
    "Avg length": "Средняя длина",
    "Base A": "Основание A",
    "Base C": "Основание C",
    "Base G": "Основание G",
    "Base T": "Основание T",
    "Base N": "Основание N",
    "GC content": "GC-состав",
    "% with adapter": "% с адаптером",
    # QC observation labels (config/qc_rules)
    "Minimum lower quartile": "Минимальный нижний квартиль",
    "Minimum median": "Минимальная медиана",
    "Modal mean read quality": "Модальное среднее качество рида",
    "Maximum A/T or G/C difference": "Макс. разница A/T или G/C",
    "Deviation from modeled distribution": "Отклонение от модельного распределения",
    "Maximum N content": "Макс. доля N",
    "Observed read lengths": "Число различных длин ридов",
    "Minimum read length": "Минимальная длина рида",
    "Sequences removed by deduplication": "Удаляется при дедупликации",
    "Maximum adapter content": "Макс. доля адаптеров",
    "All evaluated observations are within configured thresholds.": "Все наблюдения в пределах заданных порогов.",
}

CHART = {
    "Position in read (bp)": "Позиция в риде (п.н.)",
    "Phred quality": "Качество Phred",
    "Mean": "Среднее",
    "Median": "Медиана",
    "Lower quartile (Q1)": "Нижний квартиль (Q1)",
    "Mean Phred quality per read": "Среднее качество Phred рида",
    "Reads (%)": "Риды (%)",
    "NeoQC (rounded)": "NeoQC (округлено)",
    "Read distribution": "Распределение ридов",
    "FastQC-compatible (truncated)": "Совместимо с FastQC (усечено)",
    "Reads with adapter (%)": "Риды с адаптером (%)",
    "No adapters detected": "Адаптеры не обнаружены",
    "Base content (%)": "Доля оснований (%)",
    "GC content (%)": "GC-состав (%)",
    "Reads": "Риды",
    "Observed (FastQC-compatible)": "Наблюдаемое (как в FastQC)",
    "Raw observed": "Наблюдаемое (исходное)",
    "Theoretical distribution": "Теоретическое распределение",
    "N content (%)": "Доля N (%)",
    "Read length (bp)": "Длина рида (п.н.)",
    "Sequence duplication level": "Уровень дупликации",
    "Sequences (%)": "Последовательности (%)",
    "Total sequences": "Все последовательности",
    "Deduplicated sequences": "После дедупликации",
}
CHART_PATTERNS = [
    (re.compile(r"^Mean GC: (.+)$"), r"Среднее GC: \1"),
    (re.compile(r"^Mean Q(\S+)  •  Median Q(\S+)  •  Mode Q(\S+)$"), r"Среднее Q\1  •  Медиана Q\2  •  Мода Q\3"),
    (re.compile(r"^Mean (\S+) bp  •  Median (\S+) bp  •  Mode (\S+) bp$"),
     r"Среднее \1 п.н.  •  Медиана \2 п.н.  •  Мода \3 п.н."),
    (re.compile(r"^Fixed length  •  (\S+) bp$"), r"Фиксированная длина  •  \1 п.н."),
    (re.compile(r"^Peak (\S+)%$"), r"Пик \1%"),
    (re.compile(r"^(\S+) reads$"), r"\1 ридов"),
]

_REASON = re.compile(r"^(?P<label>.+?): (?P<obs>.+?); (?P<lvl>WARNING|FAIL) threshold (?P<rest>.+)\.$")


def chart_ru(text: str) -> str:
    """Translate a chart label (axis, legend or annotation) into Russian."""
    if text in CHART:
        return CHART[text]
    for pattern, repl in CHART_PATTERNS:
        if pattern.match(text):
            return pattern.sub(repl, text)
    return text


def ru(text: str) -> str:
    """Translate report UI text or a QC reason message into Russian."""
    if text in UI:
        return UI[text]
    m = _REASON.match(text)
    if m:
        label = UI.get(m["label"], m["label"])
        return f"{label}: {m['obs']}; порог {m['lvl']} {m['rest']}."
    return text
