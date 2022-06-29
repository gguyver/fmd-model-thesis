//
// Created by glen on 03/02/2020.
//

#ifndef MODEL_SINGLESEROTYPE_MATRIX_H
#define MODEL_SINGLESEROTYPE_MATRIX_H

#include <vector>
#include <fstream>

template<class T>
class Matrix {

private:
    size_t m_numRows{0}, m_numCols{0}, m_numData{0}, traversed{0};
    std::vector<T> m_Data;

    [[nodiscard]] size_t calculateIndex(size_t row, size_t col) const {
        return (row * m_numCols) + col;
    }

public:
    explicit Matrix(const T &initialVal, size_t numRows = 0, size_t numCols = 0) : m_numRows{numRows},
                                                                                   m_numCols{numCols},
                                                                                   m_numData{numRows * numCols} {
        m_Data = std::vector<T>(m_numData, initialVal);
    }

    explicit Matrix(const size_t numRows = 0, const size_t numCols = 0) : m_numRows{numRows}, m_numCols{numCols},
                                                                          m_numData{numRows * numCols},
                                                                          m_Data{std::vector<T>(numRows * numCols)} {}

    T &operator()(size_t row, size_t col) { // returns reference so can change value
        size_t index = calculateIndex(row, col);
        return m_Data.at(index);
    }

    T operator()(size_t row, size_t col) const { // const because returns a copy of the value stored
        size_t index = calculateIndex(row, col);
        return m_Data.at(index);
    }

    void addRow(const std::vector<T> &elements) {
        if (elements.size() != m_numCols) {
            throw std::invalid_argument(
                    "Matrix::addRow() must be provided a vector of the same length as the number of columns");
        }
        for (const auto &elem : elements) {
            m_Data[traversed] = elem;
            ++traversed;
        }
    }

    void writeMatrixToFile(const std::string &fileName) const {
        std::ofstream outputLocation{fileName.c_str()};
        if (outputLocation.is_open()) {
            for (size_t row = 0; row < m_numRows; ++row) {
                for (size_t col = 0; col < m_numCols; ++col) {
                    size_t index{calculateIndex(row, col)};
                    outputLocation << m_Data[index];
                    // comma/newline
                    if (col != m_numCols - 1) {
                        outputLocation << ',';
                    }
                }
                outputLocation << '\n';
            }
            outputLocation.close();
        } else {
            throw std::ofstream::failure("Cannot write Matrix to " + fileName);
        }
    }

    [[nodiscard]]  size_t getNumRows() const { return m_numRows; };

    [[nodiscard]]  size_t getNumCols() const { return m_numCols; };

    size_t getNumData() { return m_numData; };
};


#endif //MODEL_SINGLESEROTYPE_MATRIX_H
