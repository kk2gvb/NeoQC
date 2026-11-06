#include "per_base_metrics.h"

#include <algorithm>
#include <stdexcept>

namespace {

template <typename T>
void addVector(std::vector<T>& destination, const std::vector<T>& source)
{
    if (destination.size() < source.size()) {
        destination.resize(source.size(), T{});
    }
    for (std::size_t index = 0; index < source.size(); ++index) {
        destination[index] += source[index];
    }
}

} // namespace

void PerBaseMetrics::resize(std::size_t positions)
{
    baseCountA.resize(positions, 0);
    baseCountC.resize(positions, 0);
    baseCountG.resize(positions, 0);
    baseCountT.resize(positions, 0);
    baseCountN.resize(positions, 0);
    readsPerPosition.resize(positions, 0);
    qualitySum.resize(positions, 0);
    qualityCount.resize(positions, 0);
    qualityHistogram.resize(positions);
}

void PerBaseMetrics::merge(const PerBaseMetrics& other)
{
    addVector(baseCountA, other.baseCountA);
    addVector(baseCountC, other.baseCountC);
    addVector(baseCountG, other.baseCountG);
    addVector(baseCountT, other.baseCountT);
    addVector(baseCountN, other.baseCountN);
    addVector(readsPerPosition, other.readsPerPosition);
    addVector(qualitySum, other.qualitySum);
    addVector(qualityCount, other.qualityCount);

    if (qualityHistogram.size() < other.qualityHistogram.size()) {
        qualityHistogram.resize(other.qualityHistogram.size());
    }
    for (std::size_t position = 0; position < other.qualityHistogram.size(); ++position) {
        for (std::size_t quality = 0; quality < qualityHistogram[position].size(); ++quality) {
            qualityHistogram[position][quality] += other.qualityHistogram[position][quality];
        }
    }
}

PerBaseMetrics collectPerBaseMetricsCpu(
    const std::vector<FastqRecord>& records,
    const std::vector<std::uint8_t>* selected)
{
    if (selected != nullptr && selected->size() != records.size()) {
        throw std::runtime_error("Per-base metrics selection mask has an invalid size.");
    }

    PerBaseMetrics metrics;
    for (std::size_t recordIndex = 0; recordIndex < records.size(); ++recordIndex) {
        if (selected != nullptr && !(*selected)[recordIndex]) {
            continue;
        }

        const FastqRecord& record = records[recordIndex];
        const std::size_t length = record.sequence.size();
        if (length == 0) {
            continue;
        }
        metrics.resize(std::max(metrics.readsPerPosition.size(), length));

        for (std::size_t position = 0; position < length; ++position) {
            ++metrics.readsPerPosition[position];
            switch (record.sequence[position]) {
            case 'A': case 'a': ++metrics.baseCountA[position]; break;
            case 'C': case 'c': ++metrics.baseCountC[position]; break;
            case 'G': case 'g': ++metrics.baseCountG[position]; break;
            case 'T': case 't': ++metrics.baseCountT[position]; break;
            case 'N': case 'n': ++metrics.baseCountN[position]; break;
            default: break; // Matches QualityAnalyzer: invalid bases are not positional N.
            }

            if (position < record.quality.size()) {
                const int quality = static_cast<int>(record.quality[position]) - 33;
                if (quality >= 0) {
                    metrics.qualitySum[position] += static_cast<std::uint64_t>(quality);
                    ++metrics.qualityCount[position];
                    if (quality < 94) {
                        ++metrics.qualityHistogram[position][static_cast<std::size_t>(quality)];
                    }
                }
            }
        }
    }
    return metrics;
}
