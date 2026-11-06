#include "quality_analyzer.h"
#include "per_base_metrics.h"

#include <cassert>
#include <vector>

namespace {

FastqRecord makeRecord(std::string sequence, std::string quality, std::uint64_t number)
{
    return {"@READ_" + std::to_string(number), std::move(sequence), "+", std::move(quality), number};
}

void assertSamePerBaseStats(const QualityStats& expected, const QualityStats& actual)
{
    assert(expected.baseCountA == actual.baseCountA);
    assert(expected.baseCountC == actual.baseCountC);
    assert(expected.baseCountG == actual.baseCountG);
    assert(expected.baseCountT == actual.baseCountT);
    assert(expected.baseCountN == actual.baseCountN);
    assert(expected.readsPerPosition == actual.readsPerPosition);
    assert(expected.perBaseQualityGroups.size() == actual.perBaseQualityGroups.size());
    for (std::size_t index = 0; index < expected.perBaseQualityGroups.size(); ++index) {
        const auto& left = expected.perBaseQualityGroups[index];
        const auto& right = actual.perBaseQualityGroups[index];
        assert(left.start == right.start);
        assert(left.end == right.end);
        assert(left.mean == right.mean);
        assert(left.lowerQuartile == right.lowerQuartile);
        assert(left.median == right.median);
        assert(left.evaluated == right.evaluated);
    }
}

} // namespace

int main()
{
    const std::vector<FastqRecord> records = {
        makeRecord("ACGTN", "!5I~+", 1),
        makeRecord("acg", "!~I", 2),
        makeRecord("NNTA", "I!5~", 3),
    };

    QualityAnalyzer reference;
    QualityAnalyzer deferred;
    for (const FastqRecord& record : records) {
        assert(!reference.processRecord(record).found);
        assert(!deferred.processRecord(record, false).found);
    }
    deferred.applyPerBaseMetrics(collectPerBaseMetricsCpu(records));
    assertSamePerBaseStats(reference.getStats(), deferred.getStats());

    const std::vector<std::uint8_t> selected = {1, 0, 1};
    QualityAnalyzer selectedReference;
    QualityAnalyzer selectedDeferred;
    for (std::size_t index = 0; index < records.size(); ++index) {
        if (selected[index]) {
            selectedReference.processRecord(records[index]);
            selectedDeferred.processRecord(records[index], false);
        }
    }
    selectedDeferred.applyPerBaseMetrics(collectPerBaseMetricsCpu(records, &selected));
    assertSamePerBaseStats(selectedReference.getStats(), selectedDeferred.getStats());
    return 0;
}
