//
// Created by Glen on 08/05/2020.
//

#include "gtest/gtest.h"
#include "DiseaseCompartment.h"

/* Want to test that DiseaseCompartment can:
 * - Initialise properly and throw if inputs aren't correct
 * - Keep track of population properly
 * - Properly handle tau-leap
 * - Properly handle births/additions
 * - Properly handle deaths*/


// Test construction error-checking before fixtures
TEST(disease_compartment_construction, throws_on_same_statuses) {
    EXPECT_THROW(DiseaseCompartment(Status::REC, Status::REC, 100, 1), std::invalid_argument);
}

TEST(disease_compartment_construction, throws_on_invalid_subcompartments) {
    EXPECT_THROW(DiseaseCompartment(Status::REC, Status::SUS, 100, 0), std::invalid_argument);
    EXPECT_THROW(DiseaseCompartment(Status::REC, Status::SUS, 100, -1), std::invalid_argument);
}

TEST(disease_compartment_construction, throws_on_invalid_population) {
    EXPECT_THROW(DiseaseCompartment(Status::REC, Status::SUS, -1, 1), std::invalid_argument);
}

TEST(disease_compartment_construction, throws_on_invalid_compartment_switch) {
    EXPECT_THROW(DiseaseCompartment(Status::REC, CompartmentSwitch(Status::SUS, Status::REC, 1.0), 100, 1),
                 std::invalid_argument);
}

// RNG for the random functions


struct single_compartment_test : testing::Test {
    DiseaseCompartment single_compartment{Status::SUS, Status::REC, 100, 1};
};

TEST_F(single_compartment_test, population_assigned_correctly) {
    EXPECT_EQ(single_compartment.getTotalPopulation(), 100);
    EXPECT_EQ(single_compartment.getSubPopulation(0), 100);
}

TEST_F(single_compartment_test, range_check_subcompartment_access) {
    EXPECT_THROW(single_compartment.getSubPopulation(-1), std::out_of_range);
    EXPECT_THROW(single_compartment.getSubPopulation(2), std::out_of_range);
}

TEST_F(single_compartment_test, getters_work) {
    EXPECT_EQ(single_compartment.getName(), Status::SUS);
    EXPECT_EQ(std::get<Status>(single_compartment.getNextCompartment()), Status::REC);
    EXPECT_EQ(single_compartment.getNumSubcompartments(), 1);
}

TEST_F(single_compartment_test, adding_works) {
    single_compartment.addPopulation(50);
    EXPECT_EQ(single_compartment.getTotalPopulation(), 150);
    EXPECT_EQ(single_compartment.getSubPopulation(0), 150);
}

TEST_F(single_compartment_test, random_deaths_value_check) {
    std::mt19937 testPRNG(1);
    EXPECT_THROW(single_compartment.randomDeaths(testPRNG, -1), std::invalid_argument);
}

TEST_F(single_compartment_test, random_deaths_work) {
    // with a single subcompartment there is no difference between random and ordered deaths
    // also with single subcompartment, which compartments are decremented is not random
    std::mt19937 testPRNG(1);
    single_compartment.randomDeaths(testPRNG, 20);
    EXPECT_EQ(single_compartment.getTotalPopulation(), 80);
    EXPECT_EQ(single_compartment.getSubPopulation(0), 80);
}

TEST_F(single_compartment_test, random_deaths_returns_excess) {
    std::mt19937 testPRNG(1);
    EXPECT_EQ(single_compartment.randomDeaths(testPRNG, 120), 20);
}

TEST_F(single_compartment_test, ordered_deaths_value_check) {
    EXPECT_THROW(single_compartment.orderedDeaths(-1), std::invalid_argument);
}

TEST_F(single_compartment_test, ordered_deaths_work) {
    // with a single subcompartment there is no difference between random and ordered deaths
    single_compartment.orderedDeaths(20);
    EXPECT_EQ(single_compartment.getTotalPopulation(), 80);
    EXPECT_EQ(single_compartment.getSubPopulation(0), 80);
}

TEST_F(single_compartment_test, ordered_deaths_return_excess) {
    EXPECT_EQ(single_compartment.orderedDeaths(120), 20);
}

TEST_F(single_compartment_test, tau_leap_range_checks) {
    std::mt19937 testPRNG(1);
    EXPECT_THROW(single_compartment.tauLeap(testPRNG, -0.5), std::invalid_argument);
}

TEST_F(single_compartment_test, tau_leap_return_correct_value) {
    std::mt19937 testPRNG(1);
    // known fixed seed for PRNG, should always come back with 47
    auto res{single_compartment.tauLeap(testPRNG, 0.5)};
    EXPECT_EQ(std::get<Status>(res.targetCompartments), Status::REC);
    EXPECT_EQ(res.numTransferred, 47);
}

TEST_F(single_compartment_test, tau_leap_correctly_adjust_population) {
    std::mt19937 testPRNG(1);
    single_compartment.tauLeap(testPRNG, 0.5);
    EXPECT_EQ(single_compartment.getTotalPopulation(), 53);
}

struct two_compartment_test : testing::Test {
    DiseaseCompartment two_compartment{Status::SUS, Status::REC, 100, 2};
};

TEST_F(two_compartment_test, population_assigned_correctly) {
    EXPECT_EQ(two_compartment.getTotalPopulation(), 100);
    EXPECT_EQ(two_compartment.getSubPopulation(0), 100);
    EXPECT_EQ(two_compartment.getSubPopulation(1), 0);
}

TEST_F(two_compartment_test, range_check_subcompartment_access) {
    EXPECT_THROW(two_compartment.getSubPopulation(-1), std::out_of_range);
    EXPECT_NO_THROW(two_compartment.getSubPopulation(1));
    EXPECT_THROW(two_compartment.getSubPopulation(2), std::out_of_range);
}


TEST_F(two_compartment_test, adding_works) {
    two_compartment.addPopulation(50);
    EXPECT_EQ(two_compartment.getTotalPopulation(), 150);
    EXPECT_EQ(two_compartment.getSubPopulation(0), 150);
    EXPECT_EQ(two_compartment.getSubPopulation(1), 0);
}

TEST_F(two_compartment_test, random_deaths_value_check) {
    std::mt19937 testPRNG(1);
    EXPECT_THROW(two_compartment.randomDeaths(testPRNG, -1), std::invalid_argument);
}

TEST_F(two_compartment_test, random_deaths_work) {
    // bit different with > 1 subcompartment, need to put pop in second subcompartment
    std::mt19937 testPRNG(1);
    two_compartment.tauLeap(testPRNG, 0.5);
    two_compartment.randomDeaths(testPRNG, 20);
    EXPECT_EQ(two_compartment.getTotalPopulation(), 80);
    EXPECT_EQ(two_compartment.getSubPopulation(0), 4);
    EXPECT_EQ(two_compartment.getSubPopulation(1), 76);
}

TEST_F(two_compartment_test, random_deaths_returns_excess) {
    std::mt19937 testPRNG(1);
    two_compartment.tauLeap(testPRNG, 0.5);
    EXPECT_EQ(two_compartment.randomDeaths(testPRNG, 120), 20);
}

TEST_F(two_compartment_test, ordered_deaths_value_check) {
    EXPECT_THROW(two_compartment.orderedDeaths(-1), std::invalid_argument);
}

TEST_F(two_compartment_test, ordered_deaths_work) {
    // "kills" those in latest subcompartment first
    std::mt19937 testPRNG(1);
    two_compartment.tauLeap(testPRNG, 0.5);
    two_compartment.orderedDeaths(20);
    EXPECT_EQ(two_compartment.getTotalPopulation(), 80);
    EXPECT_EQ(two_compartment.getSubPopulation(0), 4);
    EXPECT_EQ(two_compartment.getSubPopulation(1), 76);
}

TEST_F(two_compartment_test, ordered_deaths_return_excess) {
    std::mt19937 testPRNG(1);
    two_compartment.tauLeap(testPRNG, 0.5);
    EXPECT_EQ(two_compartment.orderedDeaths(120), 20);
}

TEST_F(two_compartment_test, tau_leap_returns_correct_value) {
    // known fixed seed for PRNG, should always come back with 47
    std::mt19937 testPRNG(1);
    auto res1{two_compartment.tauLeap(testPRNG, 0.5)};
    auto res2{two_compartment.tauLeap(testPRNG, 0.5)};
    EXPECT_EQ(std::get<Status>(res1.targetCompartments), Status::REC);
    EXPECT_EQ(res1.numTransferred, 0);
    EXPECT_EQ(std::get<Status>(res2.targetCompartments), Status::REC);
    EXPECT_EQ(res2.numTransferred, 96);
}

TEST_F(two_compartment_test, tau_leap_correctly_adjust_population) {
    std::mt19937 testPRNG(1);
    two_compartment.tauLeap(testPRNG, 0.5);
    EXPECT_EQ(two_compartment.getTotalPopulation(), 100);
    two_compartment.tauLeap(testPRNG, 0.5);
    EXPECT_EQ(two_compartment.getTotalPopulation(), 4);
}

TEST(variant_next_compartment, tau_leap_correct) {
    std::mt19937 prng(1);
    DiseaseCompartment comp = DiseaseCompartment(Status::REC, CompartmentSwitch(Status::SUS, Status::CAR, 0.5), 100, 1);
    auto res{comp.tauLeap(prng, 0.5)};
    EXPECT_EQ(std::get<CompartmentSwitch>(res.targetCompartments).compartmentA, Status::SUS);
    EXPECT_EQ(std::get<CompartmentSwitch>(res.targetCompartments).compartmentB, Status::CAR);
    EXPECT_EQ(std::get<CompartmentSwitch>(res.targetCompartments).proportionA, 0.5);
    EXPECT_EQ(res.numTransferred, 47);
}