//
// Created by glen on 30/07/2019.
//

#include <chrono>
#include "gtest/gtest.h"
#include "Node.h"
#include "Point.h"

class Config;

class Grid;

std::mt19937 mersenne(std::chrono::high_resolution_clock::now().time_since_epoch().count());
/* Fixture for testing
struct Node_uninfect_test : testing::Test {
    Node_uninfect_test() {}
    Point p{0.0, 0.0};
    Grid* g;
    auto c{std::make_shared<Config>()};
    Node testNode{"test", p, 100, g, };
};

TEST_F(Node_uninfect_test, cattleBeginUninfected) {
    EXPECT_EQ(10, testNode.getSerotypeStatus(std::string(), 'S'));
    // expect that all new cattle are uninfected
}

TEST_F(Node_uninfect_test, infectionWorks) {
    testNode.attemptToInfect(<#initializer#>);
    EXPECT_EQ(1, testNode.getSerotypeStatus(std::string(), 'E'));
    EXPECT_TRUE(testNode.isInfected());
}

TEST_F(Node_uninfect_test, uninfectWorks) {
    testNode.attemptToInfect(<#initializer#>);
    testNode.uninfect();
    EXPECT_EQ(9, testNode.getSerotypeStatus(std::string(), 'S'));
    EXPECT_EQ(0, testNode.getSerotypeStatus(std::string(), 'E'));
    EXPECT_EQ(0, testNode.getSerotypeStatus(std::string(), 'I'));
    EXPECT_EQ(1, testNode.getSerotypeStatus(std::string(), 'R'));
    EXPECT_FALSE(testNode.isInfected());
}
 */