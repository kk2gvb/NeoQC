#include "fastq_reader.h"
#include <zlib.h>
#include <stdexcept>

FastqReader::FastqReader(const std::string& filename) {
    fileHandle = gzopen(filename.c_str(), "rb");
    if (!fileHandle) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
}

FastqReader::~FastqReader() {
    if (fileHandle) gzclose((gzFile)fileHandle);
}

bool FastqReader::readNext(FastqRecord& record) {
    char buffer[4096];
    if (gzgets((gzFile)fileHandle, buffer, sizeof(buffer)) == nullptr) return false;
    record.header = buffer;

    // sequence
    gzgets((gzFile)fileHandle, buffer, sizeof(buffer));
    record.sequence = buffer;

    // +
    gzgets((gzFile)fileHandle, buffer, sizeof(buffer));

    // quality
    gzgets((gzFile)fileHandle, buffer, sizeof(buffer));
    record.quality = buffer;

    readCount++;
    return true;
}

size_t FastqReader::getReadCount() const {
    return readCount;
}
