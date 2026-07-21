#pragma once
#include <string>
#include <vector>

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

class FastqReader {
public:
    FastqReader(const std::string& filename);
    ~FastqReader();

    bool readNext(FastqRecord& record);
    size_t getReadCount() const;

private:
    void* fileHandle;  // gzFile
    size_t readCount = 0;
};
