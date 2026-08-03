#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace neoqc::report {

// Semantic status used for metrics, sections, and conclusions. The renderer
// controls its visual representation; callers never need to provide HTML.
enum class Status {
    notEvaluated,
    passed,
    warning,
    failed
};

struct MetadataItem {
    std::string label;
    std::string value;
};

struct Metric {
    std::string label;
    std::string value;
    std::string unit;
    std::string referenceRange;
    Status status = Status::notEvaluated;
};

struct Table {
    std::string title;
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;
};

struct Section {
    std::string id;
    std::string title;
    std::string summary;
    Status status = Status::notEvaluated;
    std::vector<Metric> metrics;
    std::vector<Table> tables;
};

struct Conclusion {
    std::string text;
    Status status = Status::notEvaluated;
};

// Data-only contract between the analytical pipeline and the report module.
// All text is escaped by the renderer. Do not put HTML into these fields.
struct ReportData {
    std::string title;
    std::string reportVersion;
    std::string caseId;
    std::string generatedAt;
    std::string organization;
    std::vector<MetadataItem> metadata;
    std::vector<Section> sections;
    std::vector<Conclusion> conclusions;
    std::string disclaimer;
};

struct RenderOptions {
    std::string language = "ru";
    bool includePrintButton = true;
};

// Validates the public data contract and throws std::invalid_argument when a
// required field is absent or a table is malformed.
void validateReportData(const ReportData& report);

// Produces a standalone UTF-8 HTML document with embedded CSS and no external
// assets. The returned string is suitable for saving or further transport.
std::string renderHtmlReport(
    const ReportData& report,
    const RenderOptions& options = {}
);

// Writes via a temporary file and rename so callers do not receive a partially
// written clinical report if writing fails.
void writeHtmlReport(
    const ReportData& report,
    const std::filesystem::path& outputPath,
    const RenderOptions& options = {}
);

} // namespace neoqc::report
