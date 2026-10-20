#include "trimming/trim_report.h"

#include <fstream>
#include <stdexcept>

namespace {

void validateStats(const TrimStats& stats) {
    if (stats.total_reads != stats.passed_reads + stats.discarded_reads) {
        throw std::logic_error("Trimming statistics violate read-count invariant");
    }
    if (stats.bases_before < stats.bases_after
        || stats.bases_trimmed != stats.bases_before - stats.bases_after) {
        throw std::logic_error("Trimming statistics violate base-count invariant");
    }
    if (stats.total_pairs != stats.passed_pairs + stats.discarded_pairs) {
        throw std::logic_error("Trimming statistics violate pair-count invariant");
    }
}

}  // namespace

void writeTrimmingReport(const std::filesystem::path& outputPath,
                         const TrimStats& stats) {
    validateStats(stats);

    const std::filesystem::path temporaryPath = outputPath.string() + ".tmp";
    {
        std::ofstream output(temporaryPath, std::ios::trunc);
        if (!output) {
            throw std::runtime_error(
                "Cannot write trimming report: " + temporaryPath.string());
        }

        output << "{\n"
               << "  \"reads\": {\n"
               << "    \"total\": " << stats.total_reads << ",\n"
               << "    \"passed\": " << stats.passed_reads << ",\n"
               << "    \"discarded\": " << stats.discarded_reads << "\n"
               << "  },\n"
               << "  \"bases\": {\n"
               << "    \"before\": " << stats.bases_before << ",\n"
               << "    \"after\": " << stats.bases_after << ",\n"
               << "    \"trimmed\": " << stats.bases_trimmed << "\n"
               << "  },\n"
               << "  \"trimming\": {\n"
               << "    \"quality_trimmed_reads\": " << stats.quality_trimmed_reads << ",\n"
               << "    \"adapter_trimmed_reads\": " << stats.adapter_trimmed_reads << ",\n"
               << "    \"polyG_trimmed_reads\": " << stats.polyG_trimmed_reads << ",\n"
               << "    \"polyX_trimmed_reads\": " << stats.polyX_trimmed_reads << ",\n"
               << "    \"corrected_reads\": " << stats.corrected_reads << ",\n"
               << "    \"corrected_bases\": " << stats.corrected_bases << ",\n"
               << "    \"merged_pairs\": " << stats.merged_pairs << ",\n"
               << "    \"umi_processed_reads\": " << stats.umi_processed_reads << ",\n"
               << "    \"umi_failed_reads\": " << stats.umi_failed_reads << ",\n"
               << "    \"umi_bases\": " << stats.umi_bases << "\n"
               << "  },\n"
               << "  \"filtering\": {\n"
               << "    \"too_short_reads\": " << stats.too_short_reads << "\n"
               << "  },\n"
               << "  \"output\": {\n"
               << "    \"reads\": " << stats.output_reads << ",\n"
               << "    \"mate_discarded_reads\": " << stats.mate_discarded_reads << ",\n"
               << "    \"total_pairs\": " << stats.total_pairs << ",\n"
               << "    \"passed_pairs\": " << stats.passed_pairs << ",\n"
               << "    \"discarded_pairs\": " << stats.discarded_pairs << "\n"
               << "  },\n"
               << "  \"adapter_positions\":";

        if (stats.adapter_positions.empty()) {
            output << " {}\n}\n";
        } else {
            output << " {";
            bool first = true;
            for (const auto& [position, count] : stats.adapter_positions) {
                if (!first) output << ',';
                output << "\n    \"" << position << "\": " << count;
                first = false;
            }
            output << "\n  }\n}\n";
        }
        output.flush();
        if (!output) {
            throw std::runtime_error(
                "Cannot complete trimming report: " + temporaryPath.string());
        }
    }

    std::error_code error;
    std::filesystem::rename(temporaryPath, outputPath, error);
    if (error) {
        std::filesystem::remove(temporaryPath);
        throw std::runtime_error(
            "Cannot publish trimming report: " + error.message());
    }
}
