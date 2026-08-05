"""Data-only contract shared by analysis code and the report renderer."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Sequence


class ReportValidationError(ValueError):
    """Raised when input data cannot form a consistent report."""


class Status(str, Enum):
    NOT_EVALUATED = "not_evaluated"
    PASSED = "passed"
    WARNING = "warning"
    FAILED = "failed"


def _text(value: Any, path: str, *, required: bool = False) -> str:
    if value is None:
        result = ""
    elif isinstance(value, (str, int, float, bool)):
        result = str(value)
    else:
        raise ReportValidationError(f"{path} must be text or a scalar value")
    if required and not result.strip():
        raise ReportValidationError(f"{path} is required")
    return result


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ReportValidationError(f"{path} must be an object")
    return value


def _sequence(value: Any, path: str) -> Sequence[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ReportValidationError(f"{path} must be an array")
    return value


def _status(value: Any, path: str) -> Status:
    if value is None or value == "":
        return Status.NOT_EVALUATED
    try:
        return Status(value)
    except (TypeError, ValueError) as error:
        allowed = ", ".join(status.value for status in Status)
        raise ReportValidationError(f"{path} must be one of: {allowed}") from error


@dataclass(frozen=True)
class MetadataItem:
    label: str
    value: str

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "MetadataItem":
        return cls(
            label=_text(data.get("label"), f"{path}.label", required=True),
            value=_text(data.get("value"), f"{path}.value"),
        )


@dataclass(frozen=True)
class Metric:
    label: str
    value: str
    unit: str = ""
    reference_range: str = ""
    status: Status = Status.NOT_EVALUATED

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "Metric":
        return cls(
            label=_text(data.get("label"), f"{path}.label", required=True),
            value=_text(data.get("value"), f"{path}.value"),
            unit=_text(data.get("unit"), f"{path}.unit"),
            reference_range=_text(data.get("reference_range"), f"{path}.reference_range"),
            status=_status(data.get("status"), f"{path}.status"),
        )


@dataclass(frozen=True)
class Table:
    title: str = ""
    columns: tuple[str, ...] = ()
    rows: tuple[tuple[str, ...], ...] = ()

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "Table":
        columns = tuple(
            _text(value, f"{path}.columns[{index}]", required=True)
            for index, value in enumerate(_sequence(data.get("columns"), f"{path}.columns"))
        )
        if not columns:
            raise ReportValidationError(f"{path}.columns must contain at least one column")

        rows = []
        for row_index, raw_row in enumerate(_sequence(data.get("rows"), f"{path}.rows")):
            row_path = f"{path}.rows[{row_index}]"
            row = tuple(
                _text(value, f"{row_path}[{cell_index}]")
                for cell_index, value in enumerate(_sequence(raw_row, row_path))
            )
            if len(row) != len(columns):
                raise ReportValidationError(
                    f"{row_path} has {len(row)} cells; expected {len(columns)}"
                )
            rows.append(row)
        return cls(
            title=_text(data.get("title"), f"{path}.title"),
            columns=columns,
            rows=tuple(rows),
        )


@dataclass(frozen=True)
class Section:
    section_id: str
    title: str
    summary: str = ""
    status: Status = Status.NOT_EVALUATED
    metrics: tuple[Metric, ...] = ()
    tables: tuple[Table, ...] = ()

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "Section":
        return cls(
            section_id=_text(data.get("id"), f"{path}.id", required=True),
            title=_text(data.get("title"), f"{path}.title", required=True),
            summary=_text(data.get("summary"), f"{path}.summary"),
            status=_status(data.get("status"), f"{path}.status"),
            metrics=tuple(
                Metric.from_dict(_mapping(item, f"{path}.metrics[{index}]"),
                                 f"{path}.metrics[{index}]")
                for index, item in enumerate(
                    _sequence(data.get("metrics"), f"{path}.metrics")
                )
            ),
            tables=tuple(
                Table.from_dict(_mapping(item, f"{path}.tables[{index}]"),
                                f"{path}.tables[{index}]")
                for index, item in enumerate(
                    _sequence(data.get("tables"), f"{path}.tables")
                )
            ),
        )


@dataclass(frozen=True)
class Conclusion:
    text: str
    status: Status = Status.NOT_EVALUATED

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "Conclusion":
        return cls(
            text=_text(data.get("text"), f"{path}.text", required=True),
            status=_status(data.get("status"), f"{path}.status"),
        )


@dataclass(frozen=True)
class ReportData:
    title: str
    report_version: str
    case_id: str
    generated_at: str
    organization: str = ""
    metadata: tuple[MetadataItem, ...] = field(default_factory=tuple)
    sections: tuple[Section, ...] = field(default_factory=tuple)
    conclusions: tuple[Conclusion, ...] = field(default_factory=tuple)
    disclaimer: str = ""

    @classmethod
    def from_dict(cls, raw_data: Mapping[str, Any]) -> "ReportData":
        data = _mapping(raw_data, "report")
        report = cls(
            title=_text(data.get("title"), "report.title", required=True),
            report_version=_text(
                data.get("report_version"), "report.report_version", required=True
            ),
            case_id=_text(data.get("case_id"), "report.case_id", required=True),
            generated_at=_text(
                data.get("generated_at"), "report.generated_at", required=True
            ),
            organization=_text(data.get("organization"), "report.organization"),
            metadata=tuple(
                MetadataItem.from_dict(_mapping(item, f"report.metadata[{index}]"),
                                       f"report.metadata[{index}]")
                for index, item in enumerate(
                    _sequence(data.get("metadata"), "report.metadata")
                )
            ),
            sections=tuple(
                Section.from_dict(_mapping(item, f"report.sections[{index}]"),
                                  f"report.sections[{index}]")
                for index, item in enumerate(
                    _sequence(data.get("sections"), "report.sections")
                )
            ),
            conclusions=tuple(
                Conclusion.from_dict(_mapping(item, f"report.conclusions[{index}]"),
                                     f"report.conclusions[{index}]")
                for index, item in enumerate(
                    _sequence(data.get("conclusions"), "report.conclusions")
                )
            ),
            disclaimer=_text(data.get("disclaimer"), "report.disclaimer"),
        )
        report.validate()
        return report

    def validate(self) -> None:
        """Validate objects also created directly rather than through from_dict."""
        _text(self.title, "report.title", required=True)
        _text(self.report_version, "report.report_version", required=True)
        _text(self.case_id, "report.case_id", required=True)
        _text(self.generated_at, "report.generated_at", required=True)

        section_ids: set[str] = set()
        for section_index, section in enumerate(self.sections):
            if not section.section_id.strip():
                raise ReportValidationError(
                    f"report.sections[{section_index}].id is required"
                )
            if section.section_id in section_ids:
                raise ReportValidationError(
                    f"duplicate section id: {section.section_id}"
                )
            section_ids.add(section.section_id)
            for table_index, table in enumerate(section.tables):
                if not table.columns:
                    raise ReportValidationError(
                        f"report.sections[{section_index}].tables[{table_index}] "
                        "must have columns"
                    )
                for row_index, row in enumerate(table.rows):
                    if len(row) != len(table.columns):
                        raise ReportValidationError(
                            f"report.sections[{section_index}].tables[{table_index}]"
                            f".rows[{row_index}] has an invalid cell count"
                        )
