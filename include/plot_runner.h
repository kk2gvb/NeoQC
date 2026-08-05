#pragma once
#include <string>

struct PlotOptions {
    bool includeAdapters = true;
    bool generateSvg = true;
    bool generatePng = true;
    bool strict = true;
};

struct PlotRunResult {
    bool success = false;
    std::string manifestPath;
    std::string reportPath;
};

class PlotRunner {
public:
    // Generates every available chart in one Python process. Plot failures are
    // reported but never invalidate the underlying NeoQC analysis results.
    static PlotRunResult runAll(const std::string& resultDir,
                                const std::string& plotDir,
                                const PlotOptions& options = {});
};
