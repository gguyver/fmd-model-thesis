//
// Created by Glen on 30/03/2020.
//

#include "gtest/gtest.h"
#include "Point.h"

struct point_test : testing::Test {
    Point a{0.0, 0.0};
    Point b{0.0, 0.0};
    Point c{1.0, 0.0};
    Point d{0.0, 1.0};
};

TEST_F(point_test, points_equality_Test) {
    EXPECT_TRUE(a == b);
}

TEST_F(point_test, points_not_equal) {
    EXPECT_TRUE(a != c);
}

TEST_F(point_test, point_getters) {
    EXPECT_TRUE(c.getXCoord() == 1.0);
    EXPECT_TRUE(d.getYCoord() == 1.0);
}

TEST_F(point_test, point_setters) {
    d.setXCoord(1.0);
    EXPECT_TRUE(d.getXCoord() == 1.0);
    d.setYCoord(2.0);
    EXPECT_EQ(d.getYCoord(), 2.0);
}
