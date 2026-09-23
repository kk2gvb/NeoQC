#include "trimming/adapter_trimmer.h"

#include <utility>

AdapterTrimmer::AdapterTrimmer(TrimConfig config)
    : config_(std::move(config))
{
}

const TrimConfig& AdapterTrimmer::getConfig() const noexcept {
    return config_;
}
