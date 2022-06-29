//
// Created by glen on 03/02/2020.
//

#include <iostream>
#include <chrono>
#include <algorithm>
#include "Common.h"
#include "Config.h"

Config::Config(const std::string &configFile) {
    m_configYAML = YAML::LoadFile(configFile);
    translateYAML();
    validateParameters();
    initReporting();
    m_pRNG = std::make_unique<std::mt19937>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    m_shipmentRecords = std::make_unique<std::vector<ShipmentRecord>>();
    if (m_runShipments) {
        readShipmentFile();
    }

}

void Config::translateYAML() {
    // inputs
    m_shipmentFile = GetYAMLParam<std::string>({"inputs", "shipment_file"});
    m_savedNodeInput = GetYAMLParam<bool>({"inputs", "saved_node_input"});
    if (m_savedNodeInput) {
        m_nodeFile = GetYAMLParam<std::string>({"inputs", "saved_node_file"});
    } else {
        m_nodeFile = GetYAMLParam<std::string>({"inputs", "clean_node_file"});
    }

    // outputs
    m_outputDirectory = GetYAMLParam<std::string>({"outputs", "output_directory"});
    m_uniqueID = GetYAMLParam<std::string>({"outputs", "output_id"});
    m_extraOutputs.detailedReports = GetYAMLParam<bool>({"outputs", "detailed_reports"});
    m_extraOutputs.infectionEvents = GetYAMLParam<bool>({"outputs", "infection_events"});
    m_extraOutputs.globalCompartmentSums = GetYAMLParam<bool>({"outputs", "global_compartment_sums"});
    m_extraOutputs.nodeCompartmentSums = GetYAMLParam<bool>({"outputs", "node_compartment_sums"});
    m_extraOutputs.nodeSaveState = GetYAMLParam<bool>({"outputs", "node_save_state"});
    m_verbosity = GetYAMLParam<unsigned int>({"outputs", "verbosity"});

    // general
    m_numModelRuns = GetYAMLParam<unsigned int>({"general", "num_runs"});
    m_startDay = GetYAMLParam<int>({"general", "start_day"});
    m_endDay = GetYAMLParam<int>({"general", "end_day"});
    m_numDays = m_endDay - m_startDay + 1;

    // burn in
    m_burnIn = BurnInOptions{
            GetYAMLParam<bool>({"burn_in", "implement"}),
            GetYAMLParam<int>({"burn_in", "duration"}),
            GetYAMLParam<bool>({"burn_in", "discard"})
    };
    // model options


    m_runShipments = GetYAMLParam<bool>({"model_options", "run_shipments"});
    m_decliningCarrier = GetYAMLParam<bool>({"model_options", "declining_carrier_transmission"});
    m_maternalImmunityIsActive = GetYAMLParam<bool>({"model_options", "maternal_immunity"});
    m_forceInfection = GetYAMLParam<bool>({"model_options", "implement_force_infections"});
    m_forcingRate = GetYAMLParam<double>({"model_options", "forcing_rate"});
    m_maxNodeInfectionDuration = GetYAMLParam<unsigned int>({"model_options", "max_node_infection_duration"});

    // seeding
    size_t numSeeds{m_configYAML["seeds"].size()};
    if (numSeeds > 0) {
        for (size_t seedNum = 0; seedNum < numSeeds; ++seedNum) {
            std::string seedName{"seed_" + std::to_string(seedNum + 1)};
            SeedEvent seed{
                    m_startDay + GetYAMLParam<int>({"seeds", seedName, "days_in"}),
                    GetYAMLParam<bool>({"seeds", seedName, "fixed_seed"}),
                    GetYAMLParam<int>({"seeds", seedName, "number_random_seeds"}),
                    GetYAMLParam<std::string>({"seeds", seedName, "serotype"}),
                    GetYAMLParam<std::string>({"seeds", seedName, "fixed_seed_node"})
            };
            m_seedEvents.push_back(seed);
        }
    } else {
        if (m_savedNodeInput) {
            std::cerr
                    << "\nIgnoring lack of seeds in the assumption that the saved state provided contains infection.\n";
        } else {
            throw std::invalid_argument("At least one seed must be provided.");
        }
    }

    // Population parameters
    m_dailyBirthRate = GetYAMLParam<double>({"population_parameters", "birth_rate"});
    m_dailyNaturalMortality = GetYAMLParam<double>({"population_parameters", "mortality_rate"});

    // disease & vaccine
    try {
        // single value is more common
        m_serotypes.emplace_back(GetYAMLParam<std::string>({"disease", "serotypes"}));
    } catch (...) {
        // if it fails to find a single value, try a sequence and only fail if it can't do either.
        m_serotypes = GetYAMLParam<std::vector<std::string>>({"disease", "serotypes"});
    }


    m_shipmentFomiteSpreadProbability = GetYAMLParam<double>(
            {"disease", "shipments", "fomite_transmission_probability"});
    double vaccineDurationSum{0};
    for (const auto &serotypeName : m_serotypes) {
        std::string seroPars{serotypeName + "_parameters"};
        auto vaccineDuration = GetYAMLParam<double>({"control_options", "vaccine", serotypeName, "duration"});
        ConfigSerotypeParameters currentParameters{
                GetYAMLParam<double>({"disease", seroPars, "beta"}),
                GetYAMLParam<double>({"disease", seroPars, "sigma"}),
                GetYAMLParam<double>({"disease", seroPars, "gamma"}),
                1.0 / GetYAMLParam<double>({"disease", seroPars, "immunity_duration"}),
                1.0 / GetYAMLParam<double>({"disease", seroPars, "maternal_immunity_duration"}),
                1.0 / vaccineDuration,
                GetYAMLParam<double>({"disease", seroPars, "mortality"}),
                GetYAMLParam<double>({"disease", seroPars, "transmission"}),
                GetYAMLParam<double>({"disease", seroPars, "susceptibility"}),
                GetYAMLParam<double>({"disease", seroPars, "proportion_carriers"}),
                1.0 / GetYAMLParam<double>({"disease", seroPars, "carrier_duration"}),
                GetYAMLParam<double>({"disease", seroPars, "carrier_beta"})
        };
        m_serotypeInfectionParameters.emplace(serotypeName, currentParameters);
        // check max transmission/susceptibility - take the most transmissible serotype
        if ((currentParameters.transmission * currentParameters.susceptibility) >
            (m_maxTransmission * m_maxSusceptibility)) {
            m_maxTransmission = currentParameters.transmission;
            m_maxSusceptibility = currentParameters.susceptibility;
        }
        vaccineDurationSum += vaccineDuration;
        VaccineEfficacy efficacy = {
                GetYAMLParam<double>({"control_options", "vaccine", serotypeName, "efficacy"}),
                GetYAMLParam<double>({"control_options", "vaccine", serotypeName, "efficacy_stdev"})
        };
        m_vaccine.serotypeEfficacies.emplace(serotypeName, efficacy);
    }
    m_averageVaccineDuration =
            vaccineDurationSum /
            static_cast<double>(m_serotypes.size()); // used for deciding whether to vaccinate in reactive vaccination

    // kernel
    m_a = GetYAMLParam<double>({"kernel", "a"});
    m_scale = GetYAMLParam<double>({"kernel", "scale"});
    m_shape = GetYAMLParam<double>({"kernel", "shape"});

    // control options

    m_dayEffectiveFrom = (m_startDay + GetYAMLParam<int>({"control_options", "days_in_effective_from"}));
    // detection
    m_detectionRate = GetYAMLParam<double>({"control_options", "detection", "detection_probability"});
    m_detectionDelay = GetYAMLParam<unsigned int>({"control_options", "detection", "detection_delay"});

    // vaccine - see disease section
    m_dailyVaccinationCapacity = GetYAMLParam<int>({"control_options", "daily_vaccination_capacity"});

    // vaccination policy - reactive vaccination
    m_reactiveVaccinationImplemented = GetYAMLParam<bool>({"control_options", "reactive_vaccination", "implement"});
    m_reactiveVaccinationRadius = GetYAMLParam<double>({"control_options", "reactive_vaccination", "radius"});
    m_reactiveVaccinationCoverage = GetYAMLParam<double>({"control_options", "reactive_vaccination", "coverage"});
    m_reactiveVaccinationMinInterval = GetYAMLParam<int>({"control_options", "reactive_vaccination", "min_interval"});

    // vaccination policy - mass vaccination
    m_massVaccinationImplemented = GetYAMLParam<bool>({"control_options", "mass_vaccination", "implement"});
    m_massVaccinationInterval = GetYAMLParam<unsigned int>({"control_options", "mass_vaccination", "interval"});
    m_massVaccinationCoverage = GetYAMLParam<double>({"control_options", "mass_vaccination", "coverage"});

    // control options - movement ban policy
    m_movementBansImplemented = GetYAMLParam<bool>({"control_options", "movement_ban_policy", "implement"});
    m_movementBansRadius = GetYAMLParam<double>({"control_options", "movement_ban_policy", "radius"});
    m_movementBanDuration = GetYAMLParam<unsigned int>({"control_options", "movement_ban_policy", "duration"});
    m_movementBansCompliance = GetYAMLParam<double>({"control_options", "movement_ban_policy", "compliance"});


}

void Config::validateParameters() {
    // file paths work
    if (!common::checkFilePathWorks(m_nodeFile)) {
        throw std::invalid_argument("Cannot open/read: " + m_nodeFile);
    }
    if (!common::checkFilePathWorks(m_shipmentFile)) {
        throw std::invalid_argument("Cannot open/read: " + m_shipmentFile);
    }
    // need to figure out OS-agnostic, non std::filesystem way of checking this.
//    if(!common::checkFilePathWorks(m_outputDirectory)) {
//        throw std::invalid_argument("Cannot open/read: " + m_outputDirectory);
//    }
//    std::string slash = (m_outputDirectory.back() == '/' ? "" : "/");
//    std::string detailedReport = (m_outputDirectory + slash + "detailed-reports/");
//    if(!common::checkFilePathWorks(detailedReport)) {
//        throw std::invalid_argument("Does not exist: " + detailedReport);
//    }

    // do not necessarily need seeds if reading in saved inputs.
    for (const auto &seed : m_seedEvents) {
        if (not seed.isFixed && seed.numberRandom < 1) {
            throw std::invalid_argument("There must be at least 1 seed allowed when random seeding (number_seeds)");
        }
        if (std::find(m_serotypes.begin(), m_serotypes.end(), seed.serotype) == m_serotypes.end()) {
            std::cerr << "Warning: seedSerotype " << seed.serotype <<
                      " was not found in provided serotype names. Halting execution.\n";
            std::exit(1);
        }
    }

    // m_serotypes is not empty, doesn't contain "none"
    if (m_serotypes.empty()) {
        throw std::invalid_argument("Please provide at least one serotype name");
    }
    std::vector<std::string> invalidSerotypes{"none"};
    for (const auto &inv : invalidSerotypes) {
        if (std::find(m_serotypes.begin(), m_serotypes.end(), inv) != m_serotypes.end()) {
            // i.e. invalidSerotype has matched
            throw std::invalid_argument("Invalid serotype name: " + inv);
        }
    }

    // m_numModelRuns > 0, < INT_MAX
    if (m_numModelRuns <= 0) {
        std::cerr << "Warning: entered number of model runs is 0 or less, setting to 1.\n";
        m_numModelRuns = 1;
    }

    if (m_numModelRuns > 1 && m_extraOutputs.nodeSaveState) {
        std::cerr
                << "\nWarning: Model will only save and output the state of nodes at the end of the final model simulation.\n";
    }

    // m_endDay > (m_startDay + m_burnInDuration)
    if (m_startDay > m_endDay) {
        throw std::invalid_argument("Error: startDay greater than or equal to specified endDay.");
    }
    if (m_burnIn.implement && ((m_startDay + m_burnIn.duration) >= m_endDay)) {
        throw std::invalid_argument("Error: Duration of burn-in is greater than or equal to specified endDay.");
    }
    // m_numDays > 0
    if (m_numDays <= 0) {
        throw std::invalid_argument("Error: number of days specified is 0 or less");
    }

    //m_maxNodeInfectionDuration > 0
    if(m_maxNodeInfectionDuration <= 0) {
        throw std::invalid_argument("Error: max duration of infection in a node must be a positive, non-zero integer");
    }

    // m_burnInDuration > 0 if burn_in is true
    if (m_burnIn.implement && m_burnIn.duration <= 0) {
        throw std::invalid_argument("Error: Specified duration of burn-in is 0 or less");
    }

    // m_dailyBirthRate >= 0.0
    if (m_dailyBirthRate < 0.0) {
        m_dailyBirthRate = 0.0;
    }
    // m_dailyNaturalMortality >= 0.0
    if (m_dailyNaturalMortality < 0.0) {
        m_dailyNaturalMortality = 0.0;
    }

    if (m_shipmentFomiteSpreadProbability < 0.0 || m_shipmentFomiteSpreadProbability > 100.0) {
        throw std::invalid_argument("Error: shipment fomite spread probability must be in the range [0.0, 100.0]");
    }

    // for serotype parameters: (existence of parameters checked when reading from map
    for (auto&[serotype, values] : m_serotypeInfectionParameters) {
        // beta >= 0
        if (values.beta < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_beta, using non-negative counterpart";
            values.beta = std::abs(values.beta);
        }
        // symptomaticRate >= 0
        if (values.symptomaticRate < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_gamma, using non-negative counterpart";
            values.symptomaticRate = std::abs(values.symptomaticRate);
        }
        if (values.recoveryRate < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_mu, using non-negative counterpart";
            values.recoveryRate = std::abs(values.recoveryRate);
        }
        if (values.waningNatural < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_lambda, using non-negative counterpart";
            values.waningNatural = std::abs(values.waningNatural);
        }
        if (values.maternalImmunityDecayRate < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype
                      << "_maternalImmunityDuration, using non-negative counterpart";
            values.maternalImmunityDecayRate = std::abs(values.maternalImmunityDecayRate);
        }
        if (values.waningVaccine < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype
                      << "_vaccineDuration, using non-negative counterpart";
            values.waningVaccine = std::abs(values.waningVaccine);
        }
        if (values.associatedMortalityRate < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_mortalityRate, using non-negative counterpart";
            values.associatedMortalityRate = std::abs(values.associatedMortalityRate);
        }
        if(values.transmission < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_transmission, using non-negative counterpart";
            values.transmission = std::abs(values.transmission);
        }
        if(values.susceptibility < 0.0) {
            std::cerr << "Warning: Negative value for " << serotype << "_susceptibility, using non-negative counterpart";
            values.susceptibility = std::abs(values.susceptibility);
        }

    }

    // c just random power - needs to be 2 or greater otherwise doesn't integrate to 1
    // b is a scaling factor
    // a is also a scaling factor, can be 1 if only 1 species but can't be absorbed if more than one species without
    // looking again as won't integrate to 1

    // m_a !> 1.0,
    if (std::abs(m_a - 1.0) > 0.00000001) { // comparison of a double
        m_a = 1.0;
        std::cerr << "Parameter 'a' must be 1.0, changing to 1.0 now." << '\n';
    }
    // m_scale > 0.0
    if (m_scale < 0.0) {
        m_scale = std::abs(m_scale);
        std::cerr << "Parameter 'b' must be positive, changing to positive counterpart." << '\n';
    }
    // m_shape >= 2.0
    if (m_shape < 2.0) {
        m_shape = 2.0;
        std::cerr << "Parameter 'c' must be 2.0 or higher for the kernel to function accurately. Changing now." << '\n';
    }

    // m_forcingRate >= 0 if m_forceInfection
    if (m_forceInfection && m_forcingRate < 0.0) {
        std::cerr << "Warning: Negative value for forcingRate, will use non-negative counterpart.";
        m_forcingRate = std::fabs(m_forcingRate);
    }
    // m_detectionRate >= 0.0, <= 100.0
    if(m_detectionRate < 0.0 || m_detectionRate > 100.0) {
        throw std::out_of_range("Error: Detection rate must be between 0.0 - 100.0");
    }
    // m_detectionDelay >= 0 from unsigned

    // for serotypes -> vaccine
    for(const auto& sero : m_serotypes) {
        if(m_vaccine.serotypeEfficacies.at(sero).meanEfficacy < 0.0 ||
                m_vaccine.serotypeEfficacies.at(sero).meanEfficacy > 100.0) {
            throw std::out_of_range("Error: Vaccine efficacy must be between 0.0 - 100.0 (Serotype = " + sero + ')');
        }
        if(m_vaccine.serotypeEfficacies.at(sero).stdevEfficacy < 0.0) {
            std::cerr << "Warning: The standard deviation of the efficacy must be a positive value. Using non-negative counterpart.\n";
            m_vaccine.serotypeEfficacies.at(sero).stdevEfficacy = std::abs(m_vaccine.serotypeEfficacies.at(sero).stdevEfficacy);
        }
    }

    // m_reactiveVaccitationRadius >= 0.0
    if (m_reactiveVaccinationImplemented && m_reactiveVaccinationRadius < 0.0) {
        throw std::invalid_argument("Error: reactiveVaccinationRadius is a negative value");
    }
    // m_reactiveVaccinationCoverage in [0, 100]
    if (m_reactiveVaccinationImplemented &&
        (m_reactiveVaccinationCoverage < 0.0 || m_reactiveVaccinationCoverage > 100.0)) {
        throw std::out_of_range("Error: Coverage value of reactive vaccination must be within 0-100");
    }
    if (m_dailyVaccinationCapacity < 1) {
        throw std::invalid_argument("daily_vaccination_capacity must be a positive integer.");
    }

    // m_massVaccinationInterval >= 0.0
    if (m_massVaccinationImplemented) {
        if (m_massVaccinationInterval < 1) {
            throw std::invalid_argument("Mass vaccination coverage must be a positive integer.");
        }
        if (m_massVaccinationCoverage < 0.0 || m_massVaccinationCoverage > 100.0) {
            std::cerr
                    << "\nMass vaccination coverage must be between 0 and 100. Adjusting provided value to closest limit.\n";
            m_massVaccinationCoverage = std::clamp(m_massVaccinationCoverage, 0.0, 100.0);
        }
    }

    // m_movementBansRadius >= 0.0
    if (m_movementBansImplemented && m_movementBansRadius < 0.0) {
        throw std::invalid_argument("Error: Radius of movement ban is a negative value");
    }
    // m_movementBansCompliance in [0, 100]
    if (m_movementBansImplemented && (m_movementBansCompliance < 0.0 || m_movementBansCompliance > 100.0)) {
        throw std::out_of_range("Error: Compliance value of movement ban must be within 0-100");
    }
    // m_movementBanDuration >= 0 - already unsigned int
}

void Config::initReporting() {
    m_reporting = std::make_unique<Reporting>(this);
}

void Config::writeReports() const {
    // NODE INFECTION
    if (m_verbosity >= 1) {
        std::cout << "\n\nWriting:\n- Node Infection Report.";
    }
    m_reporting->writeReport("model-report");

    // VACCINATION SUMS
    if (m_verbosity >= 1) {
        std::cout << "\n- Vaccination Report";
    }
    m_reporting->writeReport("vaccination-report");

    // MOVEMENT BANS
    if (m_verbosity >= 1) {
        std::cout << "\n- Movement Bans Report";
    }
    m_reporting->writeReport("movement-ban-stats");

    // GLOBAL COMPARTMENT SUMS
    if (m_extraOutputs.globalCompartmentSums) {
        if (m_verbosity >= 1) {
            std::cout << "\n- Global Compartment Sums.";
        }
        m_reporting->writeReport("global-compartment-sums");
    }

    // NODE COMPARTMENT SUMS
    if (m_extraOutputs.nodeCompartmentSums) {
        if (m_verbosity >= 1) {
            std::cout << "\n- Node Compartment Sums";
        }
        m_reporting->writeReport("node-compartment-sums");
    }

    // INFECTION EVENTS
    if (m_extraOutputs.infectionEvents) {
        if (m_verbosity >= 1) {
            std::cout << "\n- Infection Events.";
        }
        m_reporting->writeReport("infection-events");
    }

    // DETAILED REPORT(s)
    if (m_extraOutputs.detailedReports) {
        if (m_verbosity >= 1) {
            std::cout << "\n- Detailed Report.\n";
        }
        m_reporting->writeReport("detailed-report");
    }

}

void Config::readShipmentFile() {
    if (m_verbosity >= 1) {
        std::cout << "Reading shipment file..." << std::endl;
    }
    std::ifstream file(m_shipmentFile.c_str());
    if (file) {
        std::string line;
        while (getline(file, line)) {
            if (!std::isalpha(line[0])) {
                std::stringstream rowstream;
                rowstream << line;
                std::vector<std::string> tokens{common::split(rowstream, ',')};
                int day{std::stoi(tokens[0])};
                std::string sourceID{tokens[1]};
                std::string targetID{tokens[2]};
                int size{std::stoi(tokens[3])};
                if (day >= m_startDay && day <= m_endDay && sourceID != targetID) {
                    m_shipmentRecords->emplace_back(day, sourceID, targetID, size);
                }
            }
        }
    } else {
        throw std::runtime_error(" Cannot open file " + m_shipmentFile);
    }
    std::sort(m_shipmentRecords->begin(), m_shipmentRecords->end(),
              [&](const ShipmentRecord &a, const ShipmentRecord &b) { return a.day < b.day; });
}

