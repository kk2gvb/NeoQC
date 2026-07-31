#include "../include/sample_sheet.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace fs = std::filesystem;

namespace {

constexpr std::array<const char*, 9> kExpectedColumns = {
    "patient_id", "sample_id", "sample_role", "material", "r1", "r2",
    "platform", "library_type", "reference"
};

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(trim(field));
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    return fields;
}

[[noreturn]] void fail(std::size_t line, const std::string& reason) {
    throw std::runtime_error("Sample sheet validation error (line "
                             + std::to_string(line) + "): " + reason);
}

std::string checkedFilePath(const std::string& value, std::size_t line,
                            const std::string& column) {
    const fs::path path(value);
    std::error_code ec;
    if (!fs::exists(path, ec) || ec) fail(line, column + " file does not exist: " + value);
    if (!fs::is_regular_file(path, ec) || ec) {
        fail(line, column + " must be a regular file: " + value);
    }

    std::ifstream input(path, std::ios::binary);
    if (!input) fail(line, column + " file is not readable: " + value);

    return fs::weakly_canonical(path, ec).string();
}

} // namespace

std::vector<SampleSheetEntry> loadAndValidateSampleSheet(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open sample sheet: " + path);

    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("Sample sheet is empty: " + path);
    }

    const auto header = splitCsvLine(line);
    if (header.size() != kExpectedColumns.size()) {
        throw std::runtime_error("Sample sheet must have exactly 9 columns");
    }
    for (std::size_t index = 0; index < header.size(); ++index) {
        if (header[index] != kExpectedColumns[index]) {
            throw std::runtime_error("Unexpected column " + std::to_string(index + 1)
                                     + ": expected '" + kExpectedColumns[index]
                                     + "', got '" + header[index] + "'");
        }
    }

    std::vector<SampleSheetEntry> entries;
    std::unordered_set<std::string> sampleKeys;
    std::unordered_map<std::string, std::string> fileOwners;
    std::unordered_map<std::string, std::string> tumorReferences;
    std::unordered_map<std::string, std::string> normalReferences;

    std::size_t lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (trim(line).empty()) continue;

        const auto fields = splitCsvLine(line);
        if (fields.size() != kExpectedColumns.size()) {
            fail(lineNumber, "expected 9 comma-separated values");
        }

        SampleSheetEntry entry {
            fields[0], fields[1], fields[2], fields[3], fields[4], fields[5],
            fields[6], fields[7], fields[8], lineNumber
        };

        if (entry.patientId.empty()) fail(lineNumber, "patient_id is required");
        if (entry.sampleId.empty()) fail(lineNumber, "sample_id is required");
        if (entry.sampleRole != "tumor" && entry.sampleRole != "normal"
            && entry.sampleRole != "tumor_rna") {
            fail(lineNumber, "sample_role must be tumor, normal, or tumor_rna");
        }
        if (entry.material != "dna" && entry.material != "rna") {
            fail(lineNumber, "material must be dna or rna");
        }
        if (entry.r1.empty()) fail(lineNumber, "r1 is required");
        if (!entry.r2.empty() && entry.r1.empty()) fail(lineNumber, "r2 requires r1");
        if (entry.platform.empty()) fail(lineNumber, "platform is required");
        if (entry.libraryType.empty()) fail(lineNumber, "library_type is required");
        if (entry.reference.empty()) fail(lineNumber, "reference is required");

        if (entry.sampleRole == "tumor_rna" && entry.material != "rna") {
            fail(lineNumber, "tumor_rna must have material rna");
        }
        if ((entry.sampleRole == "tumor" || entry.sampleRole == "normal")
            && entry.material != "dna") {
            fail(lineNumber, entry.sampleRole + " must have material dna");
        }

        const std::string sampleKey = entry.patientId + "\n" + entry.sampleId;
        if (!sampleKeys.insert(sampleKey).second) {
            fail(lineNumber, "sample_id is duplicated for patient " + entry.patientId
                             + ": " + entry.sampleId);
        }

        for (const auto& [column, file] : std::array<std::pair<const char*, std::string>, 2> {
                 {{"r1", entry.r1}, {"r2", entry.r2}}}) {
            if (file.empty()) continue;
            const std::string canonicalPath = checkedFilePath(file, lineNumber, column);
            const auto [owner, inserted] = fileOwners.emplace(
                canonicalPath, entry.patientId + "/" + entry.sampleId + "/" + column
            );
            if (!inserted) {
                fail(lineNumber, "FASTQ file is already used by " + owner->second + ": " + file);
            }
        }

        if (entry.sampleRole == "tumor") tumorReferences[entry.patientId] = entry.reference;
        if (entry.sampleRole == "normal") normalReferences[entry.patientId] = entry.reference;
        entries.push_back(std::move(entry));
    }

    if (entries.empty()) throw std::runtime_error("Sample sheet contains no samples: " + path);

    for (const auto& [patientId, tumorReference] : tumorReferences) {
        const auto normal = normalReferences.find(patientId);
        if (normal != normalReferences.end() && normal->second != tumorReference) {
            throw std::runtime_error("Sample sheet validation error: patient " + patientId
                                     + " has different references for tumor (" + tumorReference
                                     + ") and normal (" + normal->second + ") DNA");
        }
    }

    return entries;
}

std::vector<std::string> validateCaseComposition(
    const std::vector<SampleSheetEntry>& entries)
{
    struct CaseRoles {
        bool hasTumorDna = false;
        bool hasNormalDna = false;
        bool hasTumorRna = false;
    };

    std::unordered_map<std::string, CaseRoles> cases;
    for (const auto& entry : entries) {
        auto& roles = cases[entry.patientId];
        if (entry.sampleRole == "tumor") roles.hasTumorDna = true;
        if (entry.sampleRole == "normal") roles.hasNormalDna = true;
        if (entry.sampleRole == "tumor_rna") roles.hasTumorRna = true;
    }

    std::vector<std::string> warnings;
    for (const auto& [patientId, roles] : cases) {
        if (roles.hasTumorDna && !roles.hasNormalDna) {
            throw std::runtime_error(
                "Case composition error: patient " + patientId
                + " has tumor DNA but no normal DNA control"
            );
        }
        if (!roles.hasTumorRna) {
            warnings.push_back(
                "Patient " + patientId
                + ": tumor RNA is not specified; mutant-gene expression cannot be confirmed."
            );
        }
    }
    return warnings;
}
