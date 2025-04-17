#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <ctime>

// Структура за съхранение на резултата
struct AlignmentResult {
    int score;
    std::string alignedSeq1;
    std::string alignedSeq2;
    int startPos1;
    int startPos2;
};

class SmithWaterman {
private:
    // Параметри за скоровете
    int matchScore;
    int mismatchPenalty;
    int gapPenalty;

    // Матрици за изчисления
    std::vector<std::vector<int>> scoreMatrix;

    // Размер на блока за паралелизация
    int blockSize;

public:
    SmithWaterman(int match = 2, int mismatch = -1, int gap = -1, int block = 16)
        : matchScore(match), mismatchPenalty(mismatch), gapPenalty(gap), blockSize(block) {
    }

    // Основна функция за алгоритъма
    AlignmentResult align(const std::string& seq1, const std::string& seq2) {
        int m = seq1.length();
        int n = seq2.length();

        // Инициализация на матрицата
        scoreMatrix.resize(m + 1, std::vector<int>(n + 1, 0));

        int maxScore = 0;
        int maxI = 0, maxJ = 0;

        // Запълване на матрицата паралелно по диагонали
        fillMatrixParallel(seq1, seq2);

        // Намиране на максимален скор
        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <= n; j++) {
                if (scoreMatrix[i][j] > maxScore) {
                    maxScore = scoreMatrix[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        // Възстановяване на подравняването
        AlignmentResult result = backtrack(seq1, seq2, maxI, maxJ);
        return result;
    }

    // Паралелно запълване на матрицата
    void fillMatrixParallel(const std::string& seq1, const std::string& seq2) {
        int m = seq1.length();
        int n = seq2.length();

        // Изчисление на броя диагонали
        int numDiags = m + n - 1;

        // За всяка диагонал
        for (int d = 0; d < numDiags; d++) {
            int startI = std::min(d, m);
            int startJ = std::max(0, d - m);
            int diagSize = std::min(startI + 1, n - startJ);

            // Паралелизиране по диагонали
#pragma omp parallel for
            for (int k = 0; k < diagSize; k++) {
                int i = startI - k;
                int j = startJ + k;

                // Изчисляване на стойността на клетката
                int match = (seq1[i - 1] == seq2[j - 1]) ? matchScore : mismatchPenalty;
                int diagonal = scoreMatrix[i - 1][j - 1] + match;
                int left = scoreMatrix[i][j - 1] + gapPenalty;
                int up = scoreMatrix[i - 1][j] + gapPenalty;

                // Smith-Waterman никога не позволява отрицателни стойности
                scoreMatrix[i][j] = std::max(0, std::max(diagonal, std::max(left, up)));
            }
        }
    }

    // Функция за проследяване на пътя назад и създаване на подравняване
    AlignmentResult backtrack(const std::string& seq1, const std::string& seq2, int endI, int endJ) {
        AlignmentResult result;
        result.score = scoreMatrix[endI][endJ];

        std::string alignedSeq1 = "";
        std::string alignedSeq2 = "";

        int i = endI;
        int j = endJ;

        // Проследяване назад докато достигнем клетка със стойност 0
        while (i > 0 && j > 0 && scoreMatrix[i][j] > 0) {
            int match = (seq1[i - 1] == seq2[j - 1]) ? matchScore : mismatchPenalty;

            // Проверка на посоката
            if (scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + match) {
                alignedSeq1 = seq1[i - 1] + alignedSeq1;
                alignedSeq2 = seq2[j - 1] + alignedSeq2;
                i--; j--;
            }
            else if (scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gapPenalty) {
                alignedSeq1 = seq1[i - 1] + alignedSeq1;
                alignedSeq2 = "-" + alignedSeq2;
                i--;
            }
            else {
                alignedSeq1 = "-" + alignedSeq1;
                alignedSeq2 = seq2[j - 1] + alignedSeq2;
                j--;
            }
        }

        result.alignedSeq1 = alignedSeq1;
        result.alignedSeq2 = alignedSeq2;
        result.startPos1 = i;
        result.startPos2 = j;

        return result;
    }

    // Функция за печатане на скор матрицата (за дебъгване)
    void printMatrix() {
        for (size_t i = 0; i < scoreMatrix.size(); i++) {
            for (size_t j = 0; j < scoreMatrix[i].size(); j++) {
                std::cout << scoreMatrix[i][j] << "\t";
            }
            std::cout << std::endl;
        }
    }
};

// Функция за запис на резултатите във файл
void writeResultToFile(const std::string& filename, const AlignmentResult& result,
    const std::string& seq1, const std::string& seq2,
    double executionTime, int numThreads) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Грешка при отваряне на файла за запис: " << filename << std::endl;
        return;
    }

    // Запис на текущото време
    std::time_t now = std::time(nullptr);
    char timeBuffer[100];
    std::strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    outFile << "=== Smith-Waterman Alignment Results ===" << std::endl;
    outFile << "Date and Time: " << timeBuffer << std::endl;
    outFile << "Number of Threads: " << numThreads << std::endl;
    outFile << "Execution Time: " << executionTime << " ms" << std::endl;
    outFile << std::endl;

    outFile << "Sequence 1 Length: " << seq1.length() << std::endl;
    outFile << "Sequence 2 Length: " << seq2.length() << std::endl;
    outFile << std::endl;

    outFile << "Alignment Score: " << result.score << std::endl;
    outFile << "Alignment Start Position (Seq1): " << result.startPos1 << std::endl;
    outFile << "Alignment Start Position (Seq2): " << result.startPos2 << std::endl;
    outFile << std::endl;

    outFile << "Aligned Sequence 1: " << result.alignedSeq1 << std::endl;
    outFile << "Aligned Sequence 2: " << result.alignedSeq2 << std::endl;

    outFile.close();
}

// Функция за четене на последователност от файл
std::string readSequenceFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    std::string sequence;

    if (!inFile.is_open()) {
        std::cerr << "Грешка при отваряне на файла: " << filename << std::endl;
        return "";
    }

    std::string line;
    bool firstLine = true;

    while (std::getline(inFile, line)) {
        // Пропускане на заглавния ред ако е FASTA формат
        if (firstLine && line[0] == '>') {
            firstLine = false;
            continue;
        }

        // Премахване на празни пространства и нови редове
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        sequence += line;
        firstLine = false;
    }

    inFile.close();
    return sequence;
}

int main(int argc, char* argv[]) {
    // Проверка на входните аргументи
    if (argc < 3) {
        std::cerr << "Използване: " << argv[0] << " <file1> <file2> [threads] [output_file]" << std::endl;
        std::cerr << "  file1, file2 - файлове с последователности за подравняване" << std::endl;
        std::cerr << "  threads - брой нишки (по подразбиране: максимален брой)" << std::endl;
        std::cerr << "  output_file - изходен файл (по подразбиране: result.txt)" << std::endl;
        return 1;
    }

    std::string file1 = argv[1];
    std::string file2 = argv[2];

    // Конфигуриране на броя нишки
    int numThreads = (argc > 3) ? std::stoi(argv[3]) : omp_get_max_threads();
    omp_set_num_threads(numThreads);

    // Определяне на изходния файл
    std::string outputFile = (argc > 4) ? argv[4] : "result.txt";

    // Четене на последователностите
    std::string seq1 = readSequenceFromFile(file1);
    std::string seq2 = readSequenceFromFile(file2);

    if (seq1.empty() || seq2.empty()) {
        std::cerr << "Грешка при четене на последователностите." << std::endl;
        return 1;
    }

    std::cout << "Последователност 1: " << seq1.substr(0, 50) << "... (" << seq1.length() << " символа)" << std::endl;
    std::cout << "Последователност 2: " << seq2.substr(0, 50) << "... (" << seq2.length() << " символа)" << std::endl;
    std::cout << "Брой нишки: " << numThreads << std::endl;

    // Създаване на инстанция на алгоритъма
    SmithWaterman sw(2, -1, -1);  // match=2, mismatch=-1, gap=-1

    // Стартиране на таймера
    auto startTime = std::chrono::high_resolution_clock::now();

    // Изпълнение на алгоритъма
    AlignmentResult result = sw.align(seq1, seq2);

    // Спиране на таймера
    auto endTime = std::chrono::high_resolution_clock::now();
    double executionTime = std::chrono::duration<double, std::milli>(endTime - startTime).count();

    // Извеждане на резултати
    std::cout << "Време за изпълнение: " << executionTime << " ms" << std::endl;
    std::cout << "Скор: " << result.score << std::endl;
    std::cout << "Запазване на резултатите в: " << outputFile << std::endl;

    // Запис на резултатите
    writeResultToFile(outputFile, result, seq1, seq2, executionTime, numThreads);

    return 0;
}