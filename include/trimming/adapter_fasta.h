#pragma once

#include "trimming/trim_config.h"

#include <filesystem>
#include <vector>

std::vector<TrimAdapter> loadAdapterFasta(const std::filesystem::path& path);
