#pragma once

#include "trimming/trim_config.h"

// Configured component for future single-end and paired-end adapter trimming.
// Adapter matching is intentionally not exposed until TRIM-014--TRIM-016.
class AdapterTrimmer {
public:
    explicit AdapterTrimmer(TrimConfig config);

    const TrimConfig& getConfig() const noexcept;

private:
    TrimConfig config_;
};
