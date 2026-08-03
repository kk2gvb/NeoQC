#include "html_report.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace neoqc::report {
namespace {

std::string escapeHtml(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char ch : value) {
        switch (ch) {
            case '&': escaped += "&amp;"; break;
            case '<': escaped += "&lt;"; break;
            case '>': escaped += "&gt;"; break;
            case '"': escaped += "&quot;"; break;
            case '\'': escaped += "&#39;"; break;
            default: escaped += ch; break;
        }
    }
    return escaped;
}

const char* statusClass(Status status) {
    switch (status) {
        case Status::passed: return "passed";
        case Status::warning: return "warning";
        case Status::failed: return "failed";
        case Status::notEvaluated: return "neutral";
    }
    return "neutral";
}

const char* statusLabel(Status status) {
    switch (status) {
        case Status::passed: return "Пройдено";
        case Status::warning: return "Требует внимания";
        case Status::failed: return "Не пройдено";
        case Status::notEvaluated: return "Не оценено";
    }
    return "Не оценено";
}

void requireNotBlank(const std::string& value, const char* fieldName) {
    if (value.find_first_not_of(" \t\r\n") == std::string::npos) {
        throw std::invalid_argument(
            std::string("HTML report field is required: ") + fieldName
        );
    }
}

void renderStatus(std::ostringstream& output, Status status) {
    output << "<span class=\"status " << statusClass(status) << "\">"
           << statusLabel(status) << "</span>";
}

} // namespace

void validateReportData(const ReportData& report) {
    requireNotBlank(report.title, "title");
    requireNotBlank(report.reportVersion, "reportVersion");
    requireNotBlank(report.caseId, "caseId");
    requireNotBlank(report.generatedAt, "generatedAt");

    for (std::size_t sectionIndex = 0; sectionIndex < report.sections.size(); ++sectionIndex) {
        const auto& section = report.sections[sectionIndex];
        requireNotBlank(section.id, "sections[].id");
        requireNotBlank(section.title, "sections[].title");
        for (const auto& table : section.tables) {
            if (table.columns.empty()) {
                throw std::invalid_argument("HTML report table must have at least one column");
            }
            for (const auto& row : table.rows) {
                if (row.size() != table.columns.size()) {
                    throw std::invalid_argument(
                        "HTML report table row size does not match column count in section "
                        + std::to_string(sectionIndex)
                    );
                }
            }
        }
    }
}

std::string renderHtmlReport(const ReportData& report, const RenderOptions& options) {
    validateReportData(report);

    std::ostringstream output;
    output << "<!doctype html>\n<html lang=\"" << escapeHtml(options.language) << "\">\n"
           << "<head>\n<meta charset=\"utf-8\">\n"
           << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
           << "<title>" << escapeHtml(report.title) << " — "
           << escapeHtml(report.caseId) << "</title>\n"
           << R"HTML(<style>
:root{--ink:#18212b;--muted:#637080;--line:#dbe2e8;--panel:#f7f9fb;--brand:#155e75;--ok:#166534;--ok-bg:#dcfce7;--warn:#92400e;--warn-bg:#fef3c7;--fail:#991b1b;--fail-bg:#fee2e2;--neutral:#475569;--neutral-bg:#e2e8f0}
*{box-sizing:border-box}body{margin:0;background:#eef2f5;color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}.page{max-width:1080px;margin:32px auto;background:#fff;box-shadow:0 8px 28px #2634421a}.header{padding:40px 48px 32px;border-top:8px solid var(--brand);border-bottom:1px solid var(--line)}h1{font-size:30px;line-height:1.2;margin:0 0 10px}.subtitle{color:var(--muted);margin:0}.meta{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:0 32px;margin-top:28px}.meta-item{display:flex;justify-content:space-between;gap:20px;padding:9px 0;border-bottom:1px solid var(--line)}.meta-label{color:var(--muted)}main{padding:12px 48px 40px}.section{padding:28px 0;border-bottom:1px solid var(--line)}.section:last-child{border-bottom:0}.section-head{display:flex;align-items:center;justify-content:space-between;gap:20px;margin-bottom:12px}h2{font-size:21px;margin:0}h3{font-size:16px;margin:22px 0 8px}.summary{color:#334155;max-width:84ch;white-space:pre-line}.status{display:inline-block;border-radius:999px;padding:4px 10px;font-size:12px;font-weight:700;white-space:nowrap}.passed{color:var(--ok);background:var(--ok-bg)}.warning{color:var(--warn);background:var(--warn-bg)}.failed{color:var(--fail);background:var(--fail-bg)}.neutral{color:var(--neutral);background:var(--neutral-bg)}.metrics{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin-top:18px}.metric{padding:16px;border:1px solid var(--line);border-radius:8px;background:var(--panel)}.metric-label{color:var(--muted);font-size:13px}.metric-value{font-size:22px;font-weight:700;margin:5px 0}.unit{font-size:13px;font-weight:400;color:var(--muted)}.reference{font-size:12px;color:var(--muted);min-height:19px}.metric .status{margin-top:9px}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;vertical-align:top;padding:10px 12px;border:1px solid var(--line)}th{background:var(--panel)}.conclusions{margin:14px 0 0;padding:0;list-style:none}.conclusions li{display:flex;gap:12px;align-items:flex-start;padding:12px 0;border-bottom:1px solid var(--line)}.conclusions li:last-child{border:0}.disclaimer{padding:22px 48px;background:var(--panel);border-top:1px solid var(--line);color:var(--muted);font-size:12px;white-space:pre-line}.print-button{position:fixed;right:24px;bottom:24px;border:0;border-radius:8px;padding:11px 16px;background:var(--brand);color:#fff;font-weight:700;cursor:pointer}@media(max-width:700px){.page{margin:0}.header,main{padding-left:22px;padding-right:22px}.meta,.metrics{grid-template-columns:1fr}.disclaimer{padding:20px 22px}.table-wrap{overflow-x:auto}}@media print{body{background:#fff}.page{margin:0;max-width:none;box-shadow:none}.print-button{display:none}.section,.metric,table{break-inside:avoid}}
</style>
</head>
<body>
)HTML";

    output << "<article class=\"page\">\n<header class=\"header\">\n"
           << "<h1>" << escapeHtml(report.title) << "</h1>\n"
           << "<p class=\"subtitle\">Версия отчёта " << escapeHtml(report.reportVersion);
    if (!report.organization.empty()) {
        output << " · " << escapeHtml(report.organization);
    }
    output << "</p>\n<div class=\"meta\">\n"
           << "<div class=\"meta-item\"><span class=\"meta-label\">Идентификатор случая</span><strong>"
           << escapeHtml(report.caseId) << "</strong></div>\n"
           << "<div class=\"meta-item\"><span class=\"meta-label\">Дата формирования</span><strong>"
           << escapeHtml(report.generatedAt) << "</strong></div>\n";
    for (const auto& item : report.metadata) {
        output << "<div class=\"meta-item\"><span class=\"meta-label\">"
               << escapeHtml(item.label) << "</span><strong>" << escapeHtml(item.value)
               << "</strong></div>\n";
    }
    output << "</div>\n</header>\n<main>\n";

    for (const auto& section : report.sections) {
        output << "<section class=\"section\" id=\"" << escapeHtml(section.id) << "\">\n"
               << "<div class=\"section-head\"><h2>" << escapeHtml(section.title)
               << "</h2>";
        renderStatus(output, section.status);
        output << "</div>\n";
        if (!section.summary.empty()) {
            output << "<p class=\"summary\">" << escapeHtml(section.summary) << "</p>\n";
        }
        if (!section.metrics.empty()) {
            output << "<div class=\"metrics\">\n";
            for (const auto& metric : section.metrics) {
                output << "<div class=\"metric\"><div class=\"metric-label\">"
                       << escapeHtml(metric.label) << "</div><div class=\"metric-value\">"
                       << escapeHtml(metric.value);
                if (!metric.unit.empty()) {
                    output << " <span class=\"unit\">" << escapeHtml(metric.unit) << "</span>";
                }
                output << "</div><div class=\"reference\">";
                if (!metric.referenceRange.empty()) {
                    output << "Референс: " << escapeHtml(metric.referenceRange);
                }
                output << "</div>";
                renderStatus(output, metric.status);
                output << "</div>\n";
            }
            output << "</div>\n";
        }
        for (const auto& table : section.tables) {
            if (!table.title.empty()) output << "<h3>" << escapeHtml(table.title) << "</h3>\n";
            output << "<div class=\"table-wrap\"><table><thead><tr>";
            for (const auto& column : table.columns) {
                output << "<th scope=\"col\">" << escapeHtml(column) << "</th>";
            }
            output << "</tr></thead><tbody>\n";
            for (const auto& row : table.rows) {
                output << "<tr>";
                for (const auto& cell : row) output << "<td>" << escapeHtml(cell) << "</td>";
                output << "</tr>\n";
            }
            output << "</tbody></table></div>\n";
        }
        output << "</section>\n";
    }

    if (!report.conclusions.empty()) {
        output << "<section class=\"section\"><div class=\"section-head\"><h2>Итоговые выводы</h2></div>"
               << "<ul class=\"conclusions\">\n";
        for (const auto& conclusion : report.conclusions) {
            output << "<li>";
            renderStatus(output, conclusion.status);
            output << "<span>" << escapeHtml(conclusion.text) << "</span></li>\n";
        }
        output << "</ul></section>\n";
    }
    output << "</main>\n";
    if (!report.disclaimer.empty()) {
        output << "<footer class=\"disclaimer\">" << escapeHtml(report.disclaimer) << "</footer>\n";
    }
    output << "</article>\n";
    if (options.includePrintButton) {
        output << "<button class=\"print-button\" type=\"button\" onclick=\"window.print()\">Печать / PDF</button>\n";
    }
    output << "</body>\n</html>\n";
    return output.str();
}

void writeHtmlReport(const ReportData& report,
                     const std::filesystem::path& outputPath,
                     const RenderOptions& options) {
    const std::string html = renderHtmlReport(report, options);
    const auto parent = outputPath.parent_path();
    if (!parent.empty()) {
        std::error_code directoryError;
        std::filesystem::create_directories(parent, directoryError);
        if (directoryError) {
            throw std::runtime_error(
                "Cannot create HTML report directory: " + directoryError.message()
            );
        }
    }

    std::filesystem::path temporary = outputPath;
    temporary += ".tmp";
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Cannot write temporary HTML report: " + temporary.string());
        }
        output.write(html.data(), static_cast<std::streamsize>(html.size()));
        if (!output) {
            throw std::runtime_error("Failed while writing HTML report: " + temporary.string());
        }
    }

    std::error_code renameError;
    std::filesystem::rename(temporary, outputPath, renameError);
    if (renameError) {
        // Windows does not replace an existing destination via rename. Removing
        // it here is safe because the complete temporary report already exists.
        std::error_code removeError;
        std::filesystem::remove(outputPath, removeError);
        renameError.clear();
        std::filesystem::rename(temporary, outputPath, renameError);
    }
    if (renameError) {
        std::error_code cleanupError;
        std::filesystem::remove(temporary, cleanupError);
        throw std::runtime_error("Cannot publish HTML report: " + renameError.message());
    }
}

} // namespace neoqc::report
