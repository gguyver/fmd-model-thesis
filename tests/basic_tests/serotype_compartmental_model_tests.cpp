
// Created by Glen on 12/05/2020.

#include "Config.h"
#include "gtest/gtest.h"
#include "SerotypeCompartmentalModel.h"

/* Want to test that SerotypeCompartmentalModel can:
 * - Range check inputs
 * - Allow access to summaries of the data (getters)
 * - setters
 * - isInfected etc.
 * - Execute tau-leap
 * - serotype-specific deaths
 * - Support vaccination
 * - Support generating shipments
 * - uninfect
 * - */

const std::string configLocation{"../../../inputs/test-config.yaml"};

struct scm_test : testing::Test {
    Config c{configLocation};
    SerotypeCompartmentalModel scm{&c, "test", 100};
};

//TEST(compartmental_model_initialisation, cannot_intialise_neg_pop) {
//    try {
//        Config c{configLocation};
//    } catch(const std::exception& e) {
//        std::cout << e.what();
//    }
//    //auto scm = SerotypeCompartmentalModel(&c, "test", 100);
//}