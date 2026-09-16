#pragma once

#include <cstddef>
#include <vector>

struct GCModelValue {
    int percentage;
    double weight;
};

class GCModel {
public:
    explicit GCModel(size_t readLength);

    const std::vector<GCModelValue>&
    getModelValues(size_t gcCount) const;

private:
    size_t readLength_;
    std::vector<std::vector<GCModelValue>> models_;
};