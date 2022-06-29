//
// Created by glen on 13/02/2020.
//

#include "Reporting.h"

#include "Config.h"
#include "Node.h"

std::ostream &operator<<(std::ostream &os, DayContext const &context) {
    os << context.run << ", " << context.day << ", " << context.serotype;
    return os;
}

std::ostream &operator<<(std::ostream &os, const InfectionEvent &event) {
    os << event.context << ", " << event.type << ", " << event.sourceEpi
       << ", " << event.targetEpi << ", " << (event.targetInfectedPrior ? "true\n" : "false\n");
    return os;
}

std::ostream &operator<<(std::ostream &os, const InfectedNodeSummary &ins) {
    os << ins.context << ", " << ins.nodeID << ", " << ins.daysInfected << '\n';
    return os;
}

std::ostream &operator<<(std::ostream &os, DailySerotypeSummary const &dss) {
    os << dss.context << ", true_incidence, " << dss.true_incidence << '\n';
    os << dss.context << ", true_prevalence, " << dss.true_prevalence << '\n';
    os << dss.context << ", detected_incidence, " << dss.detected_incidence << '\n';
    os << dss.context << ", detected_prevalence, " << dss.detected_prevalence << '\n';
    return os;
}

std::ostream &operator<<(std::ostream &os, DailyCompartmentSums const &dcs) {
    os << dcs.context << ", MAT, " << dcs.sums.at(0) << '\n';
    os << dcs.context << ", SUS, " << dcs.sums.at(1) << '\n';
    os << dcs.context << ", EXP, " << dcs.sums.at(2) << '\n';
    os << dcs.context << ", INF, " << dcs.sums.at(3) << '\n';
    os << dcs.context << ", REC, " << dcs.sums.at(4) << '\n';
    os << dcs.context << ", SUSVAC, " << dcs.sums.at(5) << '\n';
    os << dcs.context << ", RECVAC, " << dcs.sums.at(6) << '\n';
    os << dcs.context << ", CAR, " << dcs.sums.at(7) << '\n';
    return os;
}

// Daily Node Compartment Sums
std::ostream &operator<<(std::ostream &os, DailyNodeCompartmentSums const &dncs) {
    os << dncs.context << ", " << dncs.nodeID << ", MAT, " << dncs.sums.at(0) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", SUS, " << dncs.sums.at(1) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", EXP, " << dncs.sums.at(2) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", INF, " << dncs.sums.at(3) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", REC, " << dncs.sums.at(4) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", SUSVAC, " << dncs.sums.at(5) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", RECVAC, " << dncs.sums.at(6) << '\n';
    os << dncs.context << ", " << dncs.nodeID << ", CAR, " << dncs.sums.at(7) << '\n';

    return os;
}

// SerotypeSaveState
std::ostream &operator<<(std::ostream &os, SerotypeSaveState const &seroSave) {
    os << seroSave.serotypeID << ", " << seroSave.population;
    for (int compartmentSum : seroSave.compartmentSums) {
        os << ", " << compartmentSum;
    }
    return os;
}

// vaccination sums
std::ostream &operator<<(std::ostream &os, DailyVaccinationSum const &dvs) {
    os << dvs.run << ", " << dvs.day << ", nodes_vaccinated, " << dvs.nodesVaccinated << '\n';
    os << dvs.run << ", " << dvs.day << ", animals_vaccinated, " << dvs.animalsVaccinated << '\n';
    return os;
}

// MovementBanStats
std::ostream &operator<<(std::ostream &os, MovementBanStats const &mbs) {
    os << mbs.run << ", " << mbs.day << ", movement_bans_ordered, " << mbs.numBansOrdered << '\n';
    os << mbs.run << ", " << mbs.day << ", movement_bans_complied, " << mbs.numBansComplied << '\n';
    os << mbs.run << ", " << mbs.day << ", movement_bans_active, " << mbs.numBansActive << '\n';
    return os;
}

Reporting::Reporting(Config *configPtr) : m_config(configPtr), m_outputDirectory(configPtr->m_outputDirectory),
                                          m_uniqueID(configPtr->m_uniqueID) {
    if (m_outputDirectory.back() != '/') {
        m_outputDirectory += "/";
    }
    initModelReport();
    initCompartmentSumsReport();
}

void Reporting::initModelReport() {
    // calculate size so can reserve
    size_t entriesPerDay{(m_config->m_serotypes.size())};
    size_t numEntries{((m_config->m_numModelRuns) * (static_cast<unsigned int>(m_config->m_numDays)) * entriesPerDay)};
    if (m_config->m_burnIn.implement) {
        numEntries += ((static_cast<unsigned long long int>(m_config->m_burnIn.duration)) * entriesPerDay);
    }
    m_nodeInfectionsReport.reserve(numEntries);

    // vaccination report is 1 number per day per run, so simpler
    m_vaccinationReport.reserve((m_config->m_numModelRuns) * (static_cast<unsigned int>(m_config->m_numDays)));
}

void Reporting::initCompartmentSumsReport() {
    size_t itemsPerDay{m_config->m_serotypes.size()};
    size_t numItems{((m_config->m_numModelRuns) * (static_cast<unsigned int>(m_config->m_numDays)) * itemsPerDay) + 1};
    if (m_config->m_burnIn.implement) {
        numItems += ((static_cast<unsigned long long int>(m_config->m_burnIn.duration)) * itemsPerDay);
    }
    m_compartmentSums.reserve(numItems);
}

void Reporting::initNodeCompartmentSumsReport(size_t const numNodes) {
    size_t itemsPerDay{m_config->m_serotypes.size() * numNodes};
    size_t numItems{((m_config->m_numModelRuns) * (static_cast<unsigned int>(m_config->m_numDays)) * itemsPerDay) + 1};
    if (m_config->m_burnIn.implement) {
        numItems += ((static_cast<unsigned long long int>(m_config->m_burnIn.duration)) * itemsPerDay);
    }
    m_compartmentSums.reserve(numItems);
}

DailySerotypeSummary
Reporting::sumSerotypeStatistics(std::vector<Node> &nodes, DayContext const &context) {

    DailySerotypeSummary summary(context);
    for (auto &nodePtr : nodes) {
        if (nodePtr.isInfected() && nodePtr.getSerotypeInfectedWith() == context.serotype) {
            ++summary.true_prevalence;
            int daysSinceInf{nodePtr.getDaysSinceInfected()};
            summary.true_incidence += (daysSinceInf == 1) ? 1 : 0;
            bool detected{nodePtr.detectInfected()};
            summary.detected_prevalence += (detected) ? 1 : 0;
            summary.detected_incidence += (detected && daysSinceInf == int(m_config->m_detectionDelay) + 1) ? 1 : 0;
        }
    }
    return summary;
}

void Reporting::updateModelReport(std::vector<Node> &nodes, const std::string &run,
                                  const std::string &day) {
    // columns: run, day, serotype, value_type, value_sum
    for (const auto &serotype : m_config->m_serotypes) {
        DailySerotypeSummary summary{sumSerotypeStatistics(nodes, DayContext(run, day, serotype))};
        m_nodeInfectionsReport.push_back(summary);
    }
}



void Reporting::updateCompartmentSums(const std::vector<Node> &nodes, const std::string &runName,
                                      const std::string &dayName) {
    // initialize with needed size first time as now have number of nodes,
    // which is not available when Reporting object first initialized
    if (m_config->m_extraOutputs.nodeCompartmentSums && m_nodeCompartmentSums.empty()) {
        initNodeCompartmentSumsReport(nodes.size());
    }

    for (const auto &s : m_config->m_serotypes) {
        DayContext context(runName, dayName, s);
        DailyCompartmentSums globalCompartmentSums(context);
        for (const auto &node : nodes) {
            DailyNodeCompartmentSums nodeCompartmentSums(context, node.getID());

            // relying on being in same order
            auto nodeSums{node.getSumCompartments(context.serotype)};
            for (size_t i = 0; i < nodeSums.size(); ++i) {
                if (m_config->m_extraOutputs.globalCompartmentSums) {
                    globalCompartmentSums.sums.at(i) += nodeSums[i];
                }
                if (m_config->m_extraOutputs.nodeCompartmentSums) {
                    nodeCompartmentSums.sums.at(i) += nodeSums[i];
                }
            }

            if (m_config->m_extraOutputs.nodeCompartmentSums) {
                m_nodeCompartmentSums.push_back(nodeCompartmentSums);
            }
        }
        if (m_config->m_extraOutputs.globalCompartmentSums) {
            m_compartmentSums.push_back(globalCompartmentSums);
        }
    }
}

void Reporting::saveNodeState(const std::vector<Node> &nodes) {
    std::string outputFile{m_outputDirectory + m_uniqueID + "-node-save-state.yaml"};
    std::ofstream outfile(outputFile);
    if (outfile.is_open()) {
        YAML::Emitter out;
        int count{1};
        out << YAML::BeginMap;
        for (const auto &node : nodes) {
            out << YAML::Key << "node_" + std::to_string(count);
            out << node;
            ++count;
        }
        out << YAML::EndMap;
        outfile << out.c_str();
        outfile.close();
    } else {
        std::cerr << "\nError: Cannot write events to: " << outputFile;
    }

}

void Reporting::addInfectionEvent(const std::string &run, const std::string &day, const std::string &type,
                                  const std::string &serotype,
                                  const std::string &sourceID, const std::string &destinationID, bool prior) {
    m_infectionEvents.emplace_back(run, day, type, serotype, sourceID, destinationID, prior);
}

void Reporting::updateDetailedReport(const std::string &run, const std::string &day, const std::vector<Node> &nodes) {
    for (const auto &node : nodes) {
        if (node.isInfected()) {
            m_dailyInfectedNodes.emplace_back(run, day, node.getID(), node.getSerotypeInfectedWith(),
                                              node.getDaysSinceInfected());
        }
    }

}

void
Reporting::addVaccinationSum(const std::string &run, const std::string &day, std::pair<int, int> vaccinationStats) {
    m_vaccinationReport.emplace_back(run, day, std::get<0>(vaccinationStats), std::get<1>(vaccinationStats));
}

void Reporting::writeReport(std::string_view report) const {
    if (report == "model-report") {
        writeTidyReport<DailySerotypeSummary>("-node-infection-report.csv",
                                              "run, day, serotype, value_type, value_sum\n",
                                              m_nodeInfectionsReport);
    }
    if (report == "vaccination-report") {
        writeTidyReport<DailyVaccinationSum>("-vaccination-report.csv",
                                             "run, day, value_type, value_sum\n",
                                             m_vaccinationReport);
    }
    if (report == "movement-ban-stats") {
        writeTidyReport<MovementBanStats>("-movement-bans-report.csv",
                                          "run, day, value_type, value_sum\n",
                                          m_movementBanStatistics);
    }
    if (report == "global-compartment-sums") {
        writeTidyReport<DailyCompartmentSums>("-global-compartment-sums.csv",
                                              "run, day, serotype, compartment, value_sum\n", m_compartmentSums);
    }
    if (report == "node-compartment-sums") {
        writeTidyReport<DailyNodeCompartmentSums>("-node-compartment-sums.csv",
                                                  "run, day, serotype, node, compartment, value_sum\n",
                                                  m_nodeCompartmentSums);
    }
    if (report == "infection-events") {
        writeTidyReport<InfectionEvent>("-infection-events.csv",
                                        "run, day, serotype, context, source_node, target_node, target_infected_prior\n",
                                        m_infectionEvents);
    }
    if (report == "detailed-report") {
        writeTidyReport<InfectedNodeSummary>("-detailed-report.csv",
                                             "run, day, serotype_infected_with, node_id, days_infected\n",
                                             m_dailyInfectedNodes);
    }
}

void Reporting::addMovementBanStats(const std::vector<Node> &nodes, const std::string &run, std::string day,
                                    std::pair<int, int> banStats) {
    int totalCurrentlyBanned{0};
    for (const auto &node : nodes) {
        if (not node.getAllowedToMove()) {
            ++totalCurrentlyBanned;
        }
    }
    m_movementBanStatistics.emplace_back(run, day, std::get<0>(banStats), std::get<1>(banStats), totalCurrentlyBanned);
}








