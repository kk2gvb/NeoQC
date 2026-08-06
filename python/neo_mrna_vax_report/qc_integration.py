"""Adapt NeoQC qc_evaluation.json into the complete report data model."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Mapping, Sequence

from .models import Metric, ReportData, ReportValidationError, Section, Status, Table


STATUS_MAP = {
    "pass": Status.PASSED,
    "warning": Status.WARNING,
    "fail": Status.FAILED,
    "not_evaluated": Status.NOT_EVALUATED,
}
STATUS_LABELS = {
    "pass": "PASS",
    "warning": "WARNING",
    "fail": "FAIL",
    "not_evaluated": "NOT EVALUATED",
}
SEVERITY = {"not_evaluated": -1, "pass": 0, "warning": 1, "fail": 2}


def _mapping(value: object, path: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ReportValidationError(f"{path} must be an object")
    return value


def _sequence(value: object, path: str) -> Sequence[object]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ReportValidationError(f"{path} must be an array")
    return value


def _text(value: object, path: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ReportValidationError(f"{path} is required")
    return value.strip()


def _format_observations(raw: object, raw_checks: object, path: str) -> str:
    observations = _mapping(raw, path)
    labels: dict[str, tuple[str, str]] = {}
    for index, raw_check in enumerate(_sequence(raw_checks, f"{path}.checks")):
        check = _mapping(raw_check, f"{path}.checks[{index}]")
        name = _text(check.get("observation"), f"{path}.checks[{index}].observation")
        label = _text(check.get("label"), f"{path}.checks[{index}].label")
        unit_value = check.get("unit")
        unit = unit_value.strip() if isinstance(unit_value, str) else ""
        labels[name] = (label, unit)
    values: list[str] = []
    for name, value in observations.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ReportValidationError(f"{path}.{name} must be numeric")
        label, unit = labels.get(str(name), (str(name), ""))
        suffix = f" {unit}" if unit else ""
        values.append(f"{label}: {float(value):.4g}{suffix}")
    return "; ".join(values) or "—"


def qc_section_from_dict(
    raw: Mapping[str, object],
    *,
    section_id: str = "sequencing-quality-control",
) -> Section:
    document = _mapping(raw, "qc_evaluation")
    if document.get("schema_version") != 1:
        raise ReportValidationError("unsupported qc_evaluation schema")
    ruleset = _mapping(document.get("ruleset"), "qc_evaluation.ruleset")
    ruleset_id = _text(ruleset.get("id"), "qc_evaluation.ruleset.id")
    ruleset_version = _text(ruleset.get("version"), "qc_evaluation.ruleset.version")
    evaluations = _sequence(document.get("evaluations"), "qc_evaluation.evaluations")

    counts = {status: 0 for status in STATUS_MAP}
    rows: list[tuple[str, ...]] = []
    statuses: list[str] = []
    for index, raw_evaluation in enumerate(evaluations):
        path = f"qc_evaluation.evaluations[{index}]"
        evaluation = _mapping(raw_evaluation, path)
        metric_id = _text(evaluation.get("metric_id"), f"{path}.metric_id")
        title = _text(evaluation.get("title"), f"{path}.title")
        read = _text(evaluation.get("read"), f"{path}.read")
        status = _text(evaluation.get("qc_status"), f"{path}.qc_status")
        if status not in STATUS_MAP or read not in {"R1", "R2"}:
            raise ReportValidationError(f"{path} has an invalid status or read")
        reasons = _sequence(evaluation.get("reasons"), f"{path}.reasons")
        messages = [
            _text(_mapping(reason, f"{path}.reasons[{reason_index}]").get("message"),
                  f"{path}.reasons[{reason_index}].message")
            for reason_index, reason in enumerate(reasons)
        ]
        counts[status] += 1
        statuses.append(status)
        rows.append(
            (
                title,
                read,
                STATUS_LABELS[status],
                _format_observations(
                    evaluation.get("observations"),
                    evaluation.get("checks"),
                    f"{path}.observations",
                ),
                " ".join(messages),
                metric_id,
            )
        )

    evaluated = [status for status in statuses if status != "not_evaluated"]
    overall = max(evaluated, key=SEVERITY.__getitem__) if evaluated else "not_evaluated"
    count_metrics = tuple(
        Metric(
            label=label,
            value=str(counts[status]),
            status=STATUS_MAP[status] if counts[status] else Status.NOT_EVALUATED,
        )
        for status, label in STATUS_LABELS.items()
    )
    return Section(
        section_id=section_id,
        title="Контроль качества исходных данных",
        summary=(
            "Техническая оценка NeoQC по версионированному набору правил "
            f"{ruleset_id} v{ruleset_version}. Статусы не являются клиническим заключением."
        ),
        status=STATUS_MAP[overall],
        metrics=count_metrics,
        tables=(
            Table(
                title="Матрица проверок NeoQC",
                columns=("Модуль", "Рид", "Статус", "Наблюдения", "Пояснение", "ID"),
                rows=tuple(rows),
            ),
        ),
    )


def attach_qc_evaluation(
    report: ReportData,
    evaluation_path: str | Path,
    *,
    section_id: str = "sequencing-quality-control",
) -> ReportData:
    path = Path(evaluation_path)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ReportValidationError(f"cannot read QC evaluation {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ReportValidationError(f"invalid JSON in QC evaluation {path}: {error}") from error
    section = qc_section_from_dict(_mapping(raw, "qc_evaluation"), section_id=section_id)
    sections = list(report.sections)
    for index, existing in enumerate(sections):
        if existing.section_id == section_id:
            sections[index] = section
            break
    else:
        sections.append(section)
    integrated = replace(report, sections=tuple(sections))
    integrated.validate()
    return integrated
