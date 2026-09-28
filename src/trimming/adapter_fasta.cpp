#include "trimming/adapter_fasta.h"

#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

bool validBase(char base) {
    const char normalized = static_cast<char>(
        std::toupper(static_cast<unsigned char>(base)));
    return normalized == 'A' || normalized == 'C' || normalized == 'G'
        || normalized == 'T' || normalized == 'N';
}

void finishRecord(std::vector<TrimAdapter>& adapters, TrimAdapter& current) {
    if (current.name.empty()) {
        return;
    }
    if (current.sequence.empty()) {
        throw std::runtime_error(
            "Adapter FASTA record '" + current.name + "' has an empty sequence");
    }
    adapters.push_back(std::move(current));
    current = {};
}

}  // namespace

std::vector<TrimAdapter> loadAdapterFasta(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Cannot open adapter FASTA: " + path.string());
    }

    std::vector<TrimAdapter> adapters;
    TrimAdapter current;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            finishRecord(adapters, current);
            current.name = line.substr(1);
            if (current.name.empty()) {
                throw std::runtime_error(
                    "Adapter FASTA has an empty header at line "
                    + std::to_string(lineNumber));
            }
            continue;
        }
        if (current.name.empty()) {
            throw std::runtime_error(
                "Adapter FASTA sequence appears before a header at line "
                + std::to_string(lineNumber));
        }
        for (char base : line) {
            if (!validBase(base)) {
                throw std::runtime_error(
                    "Adapter FASTA contains an invalid base at line "
                    + std::to_string(lineNumber));
            }
            current.sequence.push_back(static_cast<char>(
                std::toupper(static_cast<unsigned char>(base))));
        }
    }
    if (input.bad()) {
        throw std::runtime_error("Cannot read adapter FASTA: " + path.string());
    }
    finishRecord(adapters, current);
    if (adapters.empty()) {
        throw std::runtime_error("Adapter FASTA contains no records: " + path.string());
    }
    return adapters;
}
