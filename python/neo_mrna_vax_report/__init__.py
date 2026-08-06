"""Public API for complete neo-mRNA-vax HTML reports."""

from .html_report import render_html_report, write_html_report
from .models import (
    Conclusion,
    MetadataItem,
    Metric,
    ReportData,
    ReportValidationError,
    Section,
    Status,
    Table,
)
from .qc_integration import attach_qc_evaluation, qc_section_from_dict

__all__ = [
    "Conclusion",
    "MetadataItem",
    "Metric",
    "ReportData",
    "ReportValidationError",
    "Section",
    "Status",
    "Table",
    "render_html_report",
    "write_html_report",
    "attach_qc_evaluation",
    "qc_section_from_dict",
]
