#pragma once

#include <string>
#include <vector>

struct AdapterConfigEntry {
    std::string name;
    std::string sequence;
};

std::vector<AdapterConfigEntry> loadAdapterConfig(
    const std::string& path
);