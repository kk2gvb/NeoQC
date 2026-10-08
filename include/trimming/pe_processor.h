#pragma once

#include "fastq_reader.h"
#include "trimming/pe_overlap.h"
#include "trimming/trim_result.h"

#include <optional>

class PeProcessor {
public:
    void correct(FastqRecord& r1,
                 FastqRecord& r2,
                 TrimResult& r1Result,
                 TrimResult& r2Result,
                 const PeOverlap& overlap) const;

    std::optional<FastqRecord> merge(const FastqRecord& r1,
                                     const FastqRecord& r2,
                                     const PeOverlap& overlap) const;
};
