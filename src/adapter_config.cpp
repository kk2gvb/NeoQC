#include "../include/adapter_config.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(const std::string& value)
{
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }

    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool isValidDNA(const std::string& sequence)
{
    if (sequence.empty()) {
        return false;
    }

    for (char base : sequence) {
        switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            break;

        default:
            return false;
        }
    }

    return true;
}

} // namespace

std::vector<AdapterConfigEntry> loadAdapterConfig(
    const std::string& path)
{
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open adapter config: " + path);
    }

    std::vector<AdapterConfigEntry> adapters;

    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;

        line = trim(line);

        // Пустые строки и комментарии.
        if (line.empty() || line.front() == '#') {
            continue;
        }

        const auto separator = line.find('\t');

        if (separator == std::string::npos) {
            throw std::runtime_error(
                "Invalid adapter config at line "
                + std::to_string(lineNumber)
                + ": expected <name>\\t<sequence>");
        }

        const std::string name =
            trim(line.substr(0, separator));

        const std::string sequence =
            trim(line.substr(separator + 1));

        if (name.empty()) {
            throw std::runtime_error(
                "Invalid adapter config at line "
                + std::to_string(lineNumber)
                + ": adapter name is empty");
        }

        if (sequence.empty()) {
            throw std::runtime_error(
                "Invalid adapter config at line "
                + std::to_string(lineNumber)
                + ": adapter sequence is empty");
        }

        if (!isValidDNA(sequence)) {
            throw std::runtime_error(
                "Invalid adapter config at line "
                + std::to_string(lineNumber)
                + ": adapter sequence contains invalid DNA symbols");
        }

        adapters.push_back({
            name,
            sequence
        });
    }

    if (adapters.empty()) {
        throw std::runtime_error(
            "Adapter config contains no adapters: " + path);
    }

    return adapters;
}