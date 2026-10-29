#include "trimming/fastq_writer.h"
#include "fastq_reader.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr std::size_t WRITE_BUFFER_SIZE = 1U << 20;

}  // namespace

FastqWriter::FastqWriter(const std::filesystem::path& path)
    : path_(path)
    , file_(gzopen(path.string().c_str(), "wb"))
{
    if (file_ == nullptr) {
        throw std::runtime_error(
            "Cannot open trimmed FASTQ for writing: " + path_.string());
    }
    buffer_.reserve(WRITE_BUFFER_SIZE);
}

FastqWriter::~FastqWriter() {
    if (file_ != nullptr) {
        if (!buffer_.empty()) {
            try {
                writeBytes(buffer_.data(), buffer_.size());
            } catch (...) {
                // Destructors cannot report write failures; close() remains
                // the API which reports them to the caller.
            }
        }
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

    const std::size_t serializedSize = record.header.size()
        + record.sequence.size() + record.separator.size() + record.quality.size() + 4;
    buffer_.reserve(buffer_.size() + serializedSize);
    buffer_.append(record.header);
    buffer_.push_back('\n');
    buffer_.append(record.sequence);
    buffer_.push_back('\n');
    buffer_.append(record.separator);
    buffer_.push_back('\n');
    buffer_.append(record.quality);
    buffer_.push_back('\n');

    if (buffer_.size() >= WRITE_BUFFER_SIZE) {
        writeBytes(buffer_.data(), buffer_.size());
        buffer_.clear();
    }
}

void FastqWriter::close() {
    if (file_ == nullptr) {
        return;
    }
    if (!buffer_.empty()) {
        writeBytes(buffer_.data(), buffer_.size());
        buffer_.clear();
    }
    gzFile handle = file_;
    file_ = nullptr;
    if (gzclose(handle) != Z_OK) {
        throw std::runtime_error(
            "Cannot finish trimmed FASTQ: " + path_.string());
    }
}
