#include "html_report.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

using namespace neoqc::report;

namespace {

ReportData makeReport() {
    ReportData report;
    report.title = "Отчёт по персонализированной мРНК-вакцине";
    report.reportVersion = "0.1";
    report.caseId = "CASE-001";
    report.generatedAt = "2026-08-03T12:00:00Z";
    report.organization = "Клиника";
    report.metadata.push_back({"Референсный геном", "GRCh38"});

    Section section;
    section.id = "sequencing-qc";
    section.title = "Контроль качества секвенирования";
    section.summary = "Все обязательные образцы прошли контроль качества.";
    section.status = Status::passed;
    section.metrics.push_back({"Доля Q30", "94.8", "%", ">= 85", Status::passed});
    section.tables.push_back({"Образцы", {"Образец", "Роль"}, {{"T-01", "Опухоль"}}});
    report.sections.push_back(std::move(section));
    report.conclusions.push_back({"Материал пригоден для анализа.", Status::passed});
    report.disclaimer = "Документ предназначен для профессионального использования.";
    return report;
}

} // namespace

int main() {
    auto report = makeReport();
    report.metadata.push_back({"Проверка экранирования", "<script>alert('x')</script> & test"});

    const std::string html = renderHtmlReport(report, {.includePrintButton = false});
    assert(html.find("<!doctype html>") != std::string::npos);
    assert(html.find("CASE-001") != std::string::npos);
    assert(html.find("Контроль качества секвенирования") != std::string::npos);
    assert(html.find("<script>alert") == std::string::npos);
    assert(html.find("&lt;script&gt;alert(&#39;x&#39;)&lt;/script&gt; &amp; test") != std::string::npos);
    assert(html.find("print-button") != std::string::npos); // CSS class remains embedded.
    assert(html.find("onclick=\"window.print()\"") == std::string::npos);

    bool invalidTableRejected = false;
    try {
        auto invalid = makeReport();
        invalid.sections.front().tables.front().rows = {{"missing second cell"}};
        renderHtmlReport(invalid);
    } catch (const std::invalid_argument&) {
        invalidTableRejected = true;
    }
    assert(invalidTableRejected);

    const auto outputPath = std::filesystem::temp_directory_path()
                            / "neoqc_html_report_test_output" / "report.html";
    writeHtmlReport(report, outputPath, {.includePrintButton = false});
    std::ifstream input(outputPath, std::ios::binary);
    const std::string written((std::istreambuf_iterator<char>(input)),
                              std::istreambuf_iterator<char>());
    assert(written == html);
    std::filesystem::remove_all(outputPath.parent_path());
}
