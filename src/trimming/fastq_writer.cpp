#include "trimming/fastq_writer.h"
#include "fastq_reader.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

FastqWriter::FastqWriter(const std::filesystem::path& path)
    : path_(path)
    , file_(gzopen(path.string().c_str(), "wb"))
{
    if (file_ == nullptr) {
        throw std::runtime_error(
            "Cannot open trimmed FASTQ for writing: " + path_.string());
    }
}

FastqWriter::~FastqWriter() {
    if (file_ != nullptr) {
        gzclose(file_);
    }
}

void FastqWriter::writeBytes(const char* data, std::size_t size) {
    while (size > 0) {
        const auto chunk = static_cast<unsigned int>(std::min<std::size_t>(
            size, static_cast<std::size_t>(std::numeric_limits<int>::max())));
        const int written = gzwrite(file_, data, chunk);
        if (written <= 0 || static_cast<unsigned int>(written) != chunk) {
            int errorNumber = Z_OK;
            const char* message = gzerror(file_, &errorNumber);
            throw std::runtime_error(
                "Cannot write trimmed FASTQ '" + path_.string() + "': "
                + (message != nullptr ? message : "unknown gzip error"));
        }
        data += chunk;
        size -= chunk;
    }
}

void FastqWriter::write(const FastqRecord& record) {
    if (file_ == nullptr) {
        throw std::logic_error("Cannot write to a closed trimmed FASTQ");
    }
    if (record.sequence.size() != record.quality.size()) {
        throw std::invalid_argument(
            "Cannot write FASTQ record with different sequence and quality lengths");
    }

    const std::string serialized = record.header + '\n'
        + record.sequence + '\n'
        + record.separator + '\n'
        + record.quality + '\n';
    writeBytes(serialized.data(), serialized.size());
}

void FastqWriter::close() {
    if (file_ == nullptr) {
        return;
    }
    gzFile handle = file_;
    file_ = nullptr;
    if (gzclose(handle) != Z_OK) {
        throw std::runtime_error(
            "Cannot finish trimmed FASTQ: " + path_.string());
    }
}
