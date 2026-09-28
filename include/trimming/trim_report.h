#pragma once

#include "trimming/trim_stats.h"

#include <filesystem>

void writeTrimmingReport(const std::filesystem::path& outputPath,
                         const TrimStats& stats);
