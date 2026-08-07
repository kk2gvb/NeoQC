#include "../include/fastq_reader.h"
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <cstring> 
#include <chrono>

FastqReader::FastqReader(const std::string& filename, bool collectTiming)
    : filename(filename), collectTiming(collectTiming) {
    fileHandle = gzopen(filename.c_str(), "rb");
    if (!fileHandle) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
}

FastqReader::~FastqReader() {
    if (fileHandle) {
        gzclose(fileHandle);
    }
}

void FastqReader::trimNewlines(std::string& s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) {
        s.pop_back();
    }
}

bool FastqReader::isValidBase(char c) {
    switch (c) {
        case 'A': case 'a':
        case 'C': case 'c':
        case 'G': case 'g':
        case 'T': case 't':
        case 'N': case 'n':
            return true;
        default:
            return false;
    }
}

bool FastqReader::readLine(std::string& line) {
    line.clear();
    const int BUFFER_SIZE = 8192;
    char buffer[BUFFER_SIZE];

    while (true) {
        if (gzgets(fileHandle, buffer, BUFFER_SIZE) == nullptr)
        {
            // Нормальный конец файла
            if (gzeof(fileHandle))
            {
                return !line.empty();
            }

            int errnum = Z_OK;
            const char* errmsg = gzerror(fileHandle, &errnum);

            throw std::runtime_error(
                "Gzip read error in " + filename + ": " + errmsg);
        }

        line += buffer;

        int errnum = 0;
        const char* errmsg = gzerror(fileHandle, &errnum);

        if (errnum != Z_OK && errnum != Z_STREAM_END)
        {
            throw std::runtime_error(
                "Gzip read error in " + filename + ": " + errmsg);
        }

        // Проверяем, есть ли в строке перевод строки
        if (line.find('\n') != std::string::npos) {
            trimNewlines(line);
            return true;
        }

        // Если буфер переполнен — продолжаем читать
        if (std::strlen(buffer) == BUFFER_SIZE - 1) {
            continue;
        }
    }
}

bool FastqReader::readNext(FastqRecord& record) {
    const auto readStart = collectTiming
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};

    // Читаем заголовок
    if (!readLine(record.header)) {
        if (collectTiming)
            timing.readAndDecompress += std::chrono::steady_clock::now() - readStart;

        // Если это самое начало файла — FASTQ пустой
        if (readCount == 0) {
            std::ostringstream oss;
            oss << "FASTQ validation error:\n"
                << "file: " << filename << "\n"
                << "reason: FASTQ file is empty";

            throw std::runtime_error(oss.str());
        }

        return false;
    }

    if (record.header.empty()) {
        if (collectTiming)
            timing.readAndDecompress += std::chrono::steady_clock::now() - readStart;

        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: blank line found between FASTQ records";

        throw std::runtime_error(oss.str());
    }

    if (collectTiming) timing.readAndDecompress += std::chrono::steady_clock::now() - readStart;
    const auto validationStart = collectTiming
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};

    // Проверяем, что заголовок начинается с @
    if (record.header[0] != '@') {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: header does not start with '@': " << record.header;
        throw std::runtime_error(oss.str());
    }

    // Читаем последовательность
    if (!readLine(record.sequence)) {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: truncated record (missing sequence)";
        throw std::runtime_error(oss.str());
    }

    if (record.sequence.empty()) {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: empty sequence";

        throw std::runtime_error(oss.str());
    }

    // Читаем разделитель (+)
    if (!readLine(record.separator)) {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: truncated record (missing separator line)";
        throw std::runtime_error(oss.str());
    }

    if (record.separator.empty() || record.separator[0] != '+') {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: separator line does not start with '+': " << record.separator;
        throw std::runtime_error(oss.str());
    }

    // Читаем качество
    if (!readLine(record.quality)) {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: truncated record (missing quality)";
        throw std::runtime_error(oss.str());
    }

    if (record.quality.empty()) {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: empty quality string";

        throw std::runtime_error(oss.str());
    }

    // Проверяем равенство длин sequence и quality
    if (record.sequence.length() != record.quality.length()) {
        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: sequence and quality lengths differ: "
            << record.sequence.length() << " != " << record.quality.length();
        throw std::runtime_error(oss.str());
    }

    // Проверяем допустимые символы в последовательности
    for (size_t i = 0; i < record.sequence.length(); ++i) {
        if (!isValidBase(record.sequence[i])) {
            std::ostringstream oss;
            oss << "FASTQ validation error:\n"
                << "file: " << filename << "\n"
                << "record: " << (readCount + 1) << "\n"
                << "reason: invalid base '" << record.sequence[i]
                << "' at position " << (i + 1);
            throw std::runtime_error(oss.str());
        }
    }

    // Нумеруем запись
    readCount++;
    record.recordNumber = readCount;

    if (collectTiming) timing.validation += std::chrono::steady_clock::now() - validationStart;

    return true;
}

bool FastqReader::readBatch(std::vector<FastqRecord>& batch,
                            std::size_t batchSize)
{
    batch.clear();

    if (batch.capacity() < batchSize)
    {
        batch.reserve(batchSize);
    };

    FastqRecord record;

    while (batch.size() < batchSize)
    {
        if (!readNext(record))
        {
            break;
        }

        batch.emplace_back(std::move(record));
        record = FastqRecord{};
    }

    return !batch.empty();
}

std::size_t FastqReader::getReadCount() const {
    return readCount;
}

const std::string& FastqReader::getFilename() const {
    return filename;
}

const FastqReaderTiming& FastqReader::getTiming() const {
    return timing;
}
