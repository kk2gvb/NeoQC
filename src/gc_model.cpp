#include "gc_model.h"

#include <algorithm>
#include <cmath>

GCModel::GCModel(size_t readLength)
    : readLength_(readLength),
      models_(readLength + 1)
{
    std::vector<int> claimingCounts(101, 0);

    // Первый проход: считаем, сколько GC-count претендуют на каждый процентный bin.
    for (size_t pos = 0; pos <= readLength_; ++pos) {
        double lowCount = static_cast<double>(pos) - 0.5;
        double highCount = static_cast<double>(pos) + 0.5;

        lowCount = std::clamp(lowCount, 0.0,
                              static_cast<double>(readLength_));
        highCount = std::clamp(highCount, 0.0,
                               static_cast<double>(readLength_));

        int lowPercentage = static_cast<int>(
            std::lround(lowCount * 100.0 / readLength_));

        int highPercentage = static_cast<int>(
            std::lround(highCount * 100.0 / readLength_));

        for (int p = lowPercentage; p <= highPercentage; ++p) {
            claimingCounts[p]++;
        }
    }

    // Второй проход: создаём веса.
    for (size_t pos = 0; pos <= readLength_; ++pos) {
        double lowCount = static_cast<double>(pos) - 0.5;
        double highCount = static_cast<double>(pos) + 0.5;

        lowCount = std::clamp(lowCount, 0.0,
                              static_cast<double>(readLength_));
        highCount = std::clamp(highCount, 0.0,
                               static_cast<double>(readLength_));

        int lowPercentage = static_cast<int>(
            std::lround(lowCount * 100.0 / readLength_));

        int highPercentage = static_cast<int>(
            std::lround(highCount * 100.0 / readLength_));

        for (int p = lowPercentage; p <= highPercentage; ++p) {
            models_[pos].push_back({
                p,
                1.0 / claimingCounts[p]
            });
        }
    }
}

const std::vector<GCModelValue>&
GCModel::getModelValues(size_t gcCount) const
{
    return models_[gcCount];
}