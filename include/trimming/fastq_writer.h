#pragma once

#include <cstddef>
#include <filesystem>

#include <zlib.h>

struct FastqRecord;

// Streaming gzip FASTQ writer. Call close() to surface compression errors.
class FastqWriter {
public:
    explicit FastqWriter(const std::filesystem::path& path);
    ~FastqWriter();

    FastqWriter(const FastqWriter&) = delete;
    FastqWriter& operator=(const FastqWriter&) = delete;

    void write(const FastqRecord& record);
    void close();

private:
    void writeBytes(const char* data, std::size_t size);

    std::filesystem::path path_;
    gzFile file_ = nullptr;
};
