#include "html_report.h"

using neoqc::report::Conclusion;
using neoqc::report::Metric;
using neoqc::report::ReportData;
using neoqc::report::Section;
using neoqc::report::Status;

int main() {
    ReportData report {
        .title = "Отчёт по персонализированной мРНК-вакцине",
        .reportVersion = "0.1",
        .caseId = "DEMO-001",
        .generatedAt = "2026-08-03T12:00:00Z",
        .organization = "Название клиники",
        .metadata = {
            {"Референсный геном", "GRCh38"},
            {"Версия пайплайна", "NeoQC 0.1"}
        },
        .sections = {
            Section {
                .id = "sequencing-qc",
                .title = "Контроль качества секвенирования",
                .summary = "Пример раздела. Значения передаются аналитическим кодом.",
                .status = Status::passed,
                .metrics = {
                    Metric {"Число прочтений", "48 210 544", "", "", Status::passed},
                    Metric {"Доля Q30", "94.8", "%", ">= 85", Status::passed},
                    Metric {"GC-состав", "48.2", "%", "ожидаемый диапазон", Status::passed}
                },
                .tables = {}
            }
        },
        .conclusions = {
            Conclusion {"Данные секвенирования пригодны для дальнейшего анализа.", Status::passed}
        },
        .disclaimer = "Демонстрационный документ. Не предназначен для принятия клинических решений."
    };

    neoqc::report::writeHtmlReport(report, "mrna_report_example.html");
}
