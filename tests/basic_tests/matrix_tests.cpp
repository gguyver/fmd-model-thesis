//
// Created by Glen on 30/03/2020.
//

#include "gtest/gtest.h"
#include "Matrix.h"


/* Area to test:
 * - That the data is stored correctly,
 * - That the data are accessed correctly,
 * - That the data are modified correctly */
struct matrix_test : testing::Test {
    Matrix<int> m{20, 3, 4};
};

TEST_F(matrix_test, matrix_dimensions_correct) {
    EXPECT_EQ(m.getNumRows(), 3);
    EXPECT_EQ(m.getNumCols(), 4);
    EXPECT_EQ(m.getNumData(), 12); // row * col
}

TEST_F(matrix_test, matrix_values_correct) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_EQ(m(i, j), 20);
        }
    }
}

TEST_F(matrix_test, matrix_values_change_correctly) {
    m(0, 2) = 9;
    EXPECT_EQ(m(0, 2), 9);
    EXPECT_NE(m(2, 0), 9);
}

TEST_F(matrix_test, matrix_throws_out_of_range) {
    EXPECT_THROW(m(3, 4), std::out_of_range);
}