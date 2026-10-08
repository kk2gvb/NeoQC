#pragma once

#include "trimming/trim_config.h"
#include "trimming/trim_result.h"

struct FastqRecord;

class UmiProcessor {
public:
    explicit UmiProcessor(TrimConfig config);

    // Returns false when the read cannot provide the configured UMI.
    bool process(FastqRecord& record, TrimResult& result) const;

private:
    TrimConfig config_;
};
