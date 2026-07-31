#pragma once

#include <cstddef>
#include <string>
#include <vector>

struct SampleSheetEntry {
    std::string patientId;
    std::string sampleId;
    std::string sampleRole;
    std::string material;
    std::string r1;
    std::string r2;
    std::string platform;
    std::string libraryType;
    std::string reference;
    std::size_t lineNumber = 0;
};

// Reads a CSV sample sheet and verifies its structure and biological metadata.
// Throws std::runtime_error when a row is invalid.
std::vector<SampleSheetEntry> loadAndValidateSampleSheet(const std::string& path);

// Checks whether each patient case has the required combination of samples.
// Returns non-blocking warnings and throws std::runtime_error for blocking errors.
std::vector<std::string> validateCaseComposition(
    const std::vector<SampleSheetEntry>& entries
);
