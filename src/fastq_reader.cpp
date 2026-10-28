#include "../include/fastq_reader.h"
#include <stdexcept>
#include <sstream>
#include <cstring>

FastqReader::FastqReader(const std::string& filename) : filename(filename) {
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

bool FastqReader::readLine(std::string& line) {
    line.clear();

    while (true) {
        if (inputPosition_ == inputEnd_) {
            const int bytesRead = gzread(
                fileHandle,
                inputBuffer_.data(),
                static_cast<unsigned int>(inputBuffer_.size()));
            if (bytesRead <= 0) {
                if (gzeof(fileHandle)) {
                    return !line.empty();
                }
                int errorNumber = Z_OK;
                const char* message = gzerror(fileHandle, &errorNumber);
                throw std::runtime_error(
                    "Gzip read error in " + filename + ": "
                    + (message != nullptr ? message : "unknown gzip error"));
            }
            inputPosition_ = 0;
            inputEnd_ = static_cast<std::size_t>(bytesRead);

            int errorNumber = Z_OK;
            const char* message = gzerror(fileHandle, &errorNumber);
            if (errorNumber != Z_OK && errorNumber != Z_STREAM_END) {
                throw std::runtime_error(
                    "Gzip read error in " + filename + ": "
                    + (message != nullptr ? message : "unknown gzip error"));
            }
        }

        const char* begin = inputBuffer_.data() + inputPosition_;
        const std::size_t available = inputEnd_ - inputPosition_;
        const auto* newline = static_cast<const char*>(
            std::memchr(begin, '\n', available));
        if (newline == nullptr) {
            line.append(begin, available);
            inputPosition_ = inputEnd_;
            continue;
        }

        line.append(begin, static_cast<std::size_t>(newline - begin));
        inputPosition_ += static_cast<std::size_t>(newline - begin) + 1;
        while (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        return true;
    }
}

bool FastqReader::readNext(FastqRecord& record) {

    // Читаем заголовок
    if (!readLine(record.header)) {

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

        std::ostringstream oss;
        oss << "FASTQ validation error:\n"
            << "file: " << filename << "\n"
            << "record: " << (readCount + 1) << "\n"
            << "reason: blank line found between FASTQ records";

        throw std::runtime_error(oss.str());
    }

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

    // Проверяем допустимый диапазон Phred+33
    for (std::size_t position = 0;
        position < record.quality.size();
        ++position)
    {
        const unsigned char qualityChar =
            static_cast<unsigned char>(record.quality[position]);

        if (qualityChar < 33 || qualityChar > 126)
        {
            std::ostringstream oss;

            oss << "FASTQ validation error:\n"
                << "file: " << filename << "\n"
                << "record: " << (readCount + 1) << "\n"
                << "position: " << (position + 1) << "\n"
                << "reason: invalid Phred+33 quality character "
                << "(ASCII " << static_cast<unsigned int>(qualityChar)
                << ")";

            throw std::runtime_error(oss.str());
        }
    }

    // Нумеруем запись
    readCount++;
    record.recordNumber = readCount;

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
