// Defining the Node class
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <map>


#include "Node.h"
#include "Grid.h"
#include "Common.h"
#include "Config.h"

Node::Node(std::string ID, Point location, int numCattle, Grid *gridPtr, const Config *config, std::mt19937 *rng)
        : m_epiID{std::move(ID)}, m_location{location}, m_config(config), m_parentGrid(gridPtr), m_pRNG(rng),
          m_numCattle(numCattle) {
    // assume always have serotypes, no serotypes caught by input checking
    for (const auto &serotype : m_config->m_serotypes) {
        m_serotypeCompartments.emplace_back(m_config, serotype, numCattle);
    }
}

Node::Node(YAML::Node const & yaml, Grid *gridPtr, const Config *config, std::mt19937 *rng)
        : m_epiID(yaml["id"].as<std::string>()),
        m_location(Point(yaml["x_coord"].as<double>(), yaml["y_coord"].as<double>())),
        m_config(config),
        m_parentGrid(gridPtr),
        m_pRNG(rng) {
    m_numCattle = yaml["num_cattle"].as<int>();
    m_isInfected = yaml["is_infected"].as<bool>();
    m_serotypeInfectedWith = yaml["serotype_infected_with"].as<std::string>();
    m_transmissionPars.first = yaml["transmission_first"].as<double>();
    m_transmissionPars.second = yaml["transmission_second"].as<double>();
    m_infectionDetectable = yaml["infection_detectable"].as<bool>();
    m_hasBeenDetected = yaml["has_been_detected"].as<bool>();
    m_daysSinceInfected = yaml["days_since_infected"].as<int>();
    m_allowedMovements = yaml["allowed_movements"].as<bool>();
    m_daysSinceMovementBan = yaml["days_since_movement_ban"].as<int>();
    m_isQueuedForVaccination = yaml["queued_for_vaccination"].as<bool>();
    m_hasBeenVaccinated = yaml["has_been_vaccinated"].as<bool>();
    m_daysSinceVaccination = yaml["days_since_vaccination"].as<int>();
    // Assume want only those serotypes defined in config.yaml, will discard any in save file that
    // do not exist in the config.yaml parameters.
    for (const auto &serotype : m_config->m_serotypes) {
        if(yaml["serotype_" + serotype]) {
            m_serotypeCompartments.emplace_back(m_config, serotype, yaml["serotype_" + serotype]);
        } else {
            throw std::runtime_error("Serotype " + serotype + " specified in config.yaml cannot be found in save file. Aborting.");
        }
    }
}

std::vector<SerotypeCompartmentalModel>::iterator Node::findSerotype(std::string_view serotype) {
    auto it{std::find_if(m_serotypeCompartments.begin(), m_serotypeCompartments.end(),
                         [&](SerotypeCompartmentalModel &scm) { return scm.getSerotype() == serotype; })};
    if (it != m_serotypeCompartments.end()) {
        return it;
    } else {
        throw std::runtime_error("Node::findSerotype() somehow cannot find specified serotype");
    }
}

std::vector<SerotypeCompartmentalModel>::const_iterator Node::findSerotype(std::string_view serotype) const {
    auto it{std::find_if(m_serotypeCompartments.begin(), m_serotypeCompartments.end(),
                         [&](const SerotypeCompartmentalModel &scm) { return scm.getSerotype() == serotype; })};
    if (it != m_serotypeCompartments.end()) {
        return it;
    } else {
        throw std::runtime_error("Node::findSerotype() somehow cannot find specified serotype");
    }
}


std::string Node::getID() const { return m_epiID; }

Point Node::getLocation() const { return m_location; }

bool Node::isInfected() const { return m_isInfected; }

void Node::update() {
    // also advance daysSinceInfected, MovementBan etc.
    checkIsInfected();
    updateDSI();
    updateDaysSinceMovementBan();
    updateDaysSinceVaccinated();
    if (m_isSentinelNode) {
        if (m_parentCell->getMaxNodeSize() != m_numCattle) {
            m_parentCell->updateMaxNodeSize(m_numCattle);
        }
    }
}


std::pair<bool, std::string> Node::attemptToInfect(std::map<std::string, double> &serotypeInfectionProbabilities) {
    // trim serotypes cannot be infected with
    for (auto it = serotypeInfectionProbabilities.begin(); it != serotypeInfectionProbabilities.end();) {
        if (canBeInfectedWithSerotype(it->first)) {
            ++it;
        } else {
            it = serotypeInfectionProbabilities.erase(it);
        }
    }
    if (!serotypeInfectionProbabilities.empty()) {
        // select between serotype based on relative probabilities of infection.
        std::string selectedSerotype{"none"};
        int numInfected;
        if (serotypeInfectionProbabilities.size() == 1) {
            selectedSerotype = serotypeInfectionProbabilities.begin()->first;
        } else {
            std::vector<double> serotypeNumbers(serotypeInfectionProbabilities.size());
            std::vector<std::string> serotypeNames(serotypeInfectionProbabilities.size());
            for (const auto&[serotype, P] : serotypeInfectionProbabilities) {
                serotypeNames.push_back(serotype);
                serotypeNumbers.push_back(P);
            }
            std::discrete_distribution<size_t> discreteDistribution(serotypeNumbers.begin(), serotypeNumbers.end());
            selectedSerotype = serotypeNames.at(discreteDistribution(*(m_pRNG)));
        }
        auto it{findSerotype(selectedSerotype)};
        /* multiply by the proportion of the population that is actually susceptible, because assuming a random mixing of the population
         * the probability that the cattle is actually susceptible needs to be taken into account. */
        double trueProb{serotypeInfectionProbabilities.at(selectedSerotype) *
                        (double((*it).getSumSusceptible()) / double((*it).getTotalPopulation()))};
        std::binomial_distribution<int> binomialDistribution((*it).getSumSusceptible(), trueProb);
        numInfected = binomialDistribution(*(m_pRNG));

        if (numInfected > 0) {
            // infect available susceptible cattle
            bool successfullyInfected = (*it).infect(numInfected);
            // only update properties if successful and Node not already infected
            if (successfullyInfected && !m_isInfected && m_serotypeInfectedWith == "none") {
                m_isInfected = true;
                m_serotypeInfectedWith = selectedSerotype;
                // decide based on detection rate whether this is detectable
                std::uniform_real_distribution<double> uniformDist(0.0, 100.0);
                double detectable = uniformDist((*(m_pRNG)));
                detectable <= m_config->m_detectionRate ? (m_infectionDetectable = true)
                                                        : (m_infectionDetectable = false);
            }
            return std::make_pair(true, selectedSerotype);
        } else {
            return std::make_pair(false, "none");
        }
    } else {
        return std::make_pair(false, "none");
    }
}


void Node::setMovement(bool newMovementPermission) {
    m_allowedMovements = newMovementPermission;
    m_daysSinceMovementBan = 0;
}

std::pair<int, int> Node::calculateBirthsAndDeaths() const {
    // Separate the events that are not associated with specific serotypes so don't do them multiple times
    // e.g. number of births has nothing to do with serotypes
    std::vector<double> rates{m_config->m_dailyBirthRate * double(m_numCattle),
                              m_config->m_dailyNaturalMortality * double(m_numCattle)};
    std::pair<int, int> events{0, 0};

    // births
    if (rates[0] > 0.0) {
        std::poisson_distribution<int> birthDist(rates[0]);
        events.first = birthDist(*(m_pRNG));
    }

    // deaths
    if (rates[1] > 0.0) {
        std::poisson_distribution<int> deathsDist(rates[1]);
        events.second = deathsDist(*(m_pRNG));
    }

    return events;
}

void Node::tauLeap() {

    // Calculate births and non-serotype deaths before resolve serotype-specific stuff, then resolve births etc.
    auto birthsAndDeaths{calculateBirthsAndDeaths()};

    // can only be 1 serotype at a time, so only 1 SCM will have infection mortality
    int serotypeDeaths{0};
    for (auto &scm : m_serotypeCompartments) {
        scm.tauLeap(*(m_pRNG));
        if (scm.getSerotype() == m_serotypeInfectedWith) {
            serotypeDeaths += scm.infectionMortality(*(m_pRNG));
        }
    }

    // harmonise serotype deaths (each SCM should have the same population at the end)
    for (auto &scm : m_serotypeCompartments) {
        if (scm.getSerotype() != m_serotypeInfectedWith) {
            if (birthsAndDeaths.second + serotypeDeaths > 0) {
                scm.imposeDeaths(*(m_pRNG), birthsAndDeaths.second + serotypeDeaths);
            }
        } else {
            if (birthsAndDeaths.second > 0) {
                scm.imposeDeaths(*(m_pRNG), birthsAndDeaths.second);
            }
        }
        if (birthsAndDeaths.first > 0) {
            scm.births(*(m_pRNG), birthsAndDeaths.first);
        }
    }
    // double check all populations are in sync
    // TODO: make this a test
    int newPop{m_serotypeCompartments.front().getTotalPopulation()};
    for (const auto &scm : m_serotypeCompartments) {
        if (scm.getTotalPopulation() != newPop) {
            throw std::runtime_error(
                    "Error, multiple SerotypeCompartmentalModels disagree on Node total population. Called from Node::tauLeap().");
        }
    }
    m_numCattle = newPop;
    updateTransmission();
}


Shipment Node::generateShipment(int numShipped) {
    // TODO: Refactor shipping to be single function interaction between Nodes
    for (const auto &scm : m_serotypeCompartments) {
        if (scm.getTotalPopulation() != m_numCattle) {
            throw std::runtime_error(
                    "Error, multiple SerotypeCompartmentalModels disagree on Node total population. Called from Node::tauLeap().");
        }
    }
    // generate the actual shipment of cattle
    Shipment shipment;
    shipment.sourceIsInfected = m_isInfected;
    shipment.infectionSerotype = m_serotypeInfectedWith;
    // ASSUMPTION: all cattle in any class are equally likely to be shipped

    for (auto &scm : m_serotypeCompartments) {
        SerotypeCompartmentalModel shipmentSCM{scm.sampleForShipment(*(m_pRNG), numShipped)};
        shipment.serotypeCompartments.push_back(shipmentSCM);
        shipment.shipmentIsInfected = shipmentSCM.getIsInfected();
    }
    // if > population, all SCM should truncate to same number and so have same population
    shipment.size = shipment.serotypeCompartments.front().getTotalPopulation();
    m_numCattle = m_serotypeCompartments.front().getTotalPopulation();

    return shipment;
}

ShipmentResult Node::acceptShipment(Shipment &incomingShipment) {
    ShipmentResult result;
    result.targetInfectedPrior = m_isInfected;

    std::uniform_real_distribution<double> unifDist(0.0, 100.0);
    result.fomiteTransmissionOccured = (incomingShipment.sourceIsInfected &&
                                        unifDist(*m_pRNG) <= m_config->m_shipmentFomiteSpreadProbability);

    if (m_isInfected) {
        if ((incomingShipment.shipmentIsInfected || result.fomiteTransmissionOccured)) {
            if (incomingShipment.infectionSerotype != m_serotypeInfectedWith) {
                for (auto &scm : incomingShipment.serotypeCompartments) {
                    if (scm.getSerotype() == incomingShipment.infectionSerotype) {
                        scm.uninfect();
                    }
                }
                incomingShipment.shipmentIsInfected = false;
                incomingShipment.sourceIsInfected = false;
                incomingShipment.infectionSerotype = "none";
            } else {
                result.targetNowInfected = true;
                // otherwise do nothing as already infected with same serotype and shipment will be added
            }
        }
    } else {
        // can now proceed with assumption that infection is valid
        if (incomingShipment.shipmentIsInfected || result.fomiteTransmissionOccured) {
            m_isInfected = true;
            result.targetNowInfected = true;
            m_serotypeInfectedWith = incomingShipment.infectionSerotype;
            double detectable = unifDist((*(m_pRNG)));
            detectable <= m_config->m_detectionRate ? (m_infectionDetectable = true)
                                                    : (m_infectionDetectable = false);
            if (result.fomiteTransmissionOccured && not incomingShipment.shipmentIsInfected) {
                // need to add infected animals as none in shipment or farm
                for (auto &scm : m_serotypeCompartments) {
                    if (scm.getSerotype() == incomingShipment.infectionSerotype && scm.getSumSusceptible() > 0) {
                        // TODO: decide on better way of deciding number infected, should be low, e^-x?
                        scm.infect(1);
                    }
                }
            }
        }
    }
    for (auto &shipmentSCM : incomingShipment.serotypeCompartments) {
        auto it{findSerotype(shipmentSCM.getSerotype())};
        (*it) += shipmentSCM;
    }
    m_numCattle += incomingShipment.size;
    return result;
}

int Node::getDaysSinceInfected() const {
    return m_daysSinceInfected;
}

void Node::updateDSI() {
    if (m_isInfected) {
        ++m_daysSinceInfected;
    } else {
        m_daysSinceInfected = 0;
    }
}

int Node::vaccinate() {
    /* ASSUMPTION: Calling this function takes place AFTER logic has been run to decide whether vaccination should take place*/
    /* ASSUMPTION: Vaccinating cattle who have maternal immunity, are already vaccinated, or are infected doesn't do anything*/

    auto vaccine = m_config->m_vaccine;
    /* Different serotypes will have different numbers vaccinated/protected, so need to take lowest as ASSUMPTION is
     * that essentially only susceptible, recovered and carrier animals vaccinated and this might be different for different serotypes.
     * Also if was set to 0, would never update with this method.
     * More realistically _some_ of the maternally immune and vaccinate cattle would be vaccinated as well, but
     * the model has no way of tracking how long individual animals are MAT or VAC so can't really know.
     * Felt it was best to track what the _model_ was actually transferring to the vaccinated compartment. */
    int vaccinesUsed{m_numCattle};
    for (auto &scm : m_serotypeCompartments) {
        // std::normal_distribution has to return double/float/long double
        std::normal_distribution<double> normDist(vaccine.serotypeEfficacies.at(scm.getSerotype()).meanEfficacy,
                                                  vaccine.serotypeEfficacies.at(scm.getSerotype()).stdevEfficacy);
        double serotypeEfficacy{normDist((*(m_pRNG))) / 100.0};
        // bounds check
        serotypeEfficacy = std::clamp(serotypeEfficacy, 0.0, 1.0);

        // Calculate proportion that would be successfully vaccinated and update numbers
        int numVaccinated = scm.vaccinate(*(m_pRNG), serotypeEfficacy);
        vaccinesUsed = std::min(vaccinesUsed, numVaccinated);
    }
    m_isQueuedForVaccination = false;
    m_hasBeenVaccinated = true;
    m_daysSinceVaccination = 0;
    return vaccinesUsed;
}

void Node::updateDaysSinceMovementBan() {
    // update days since movement ban
    // to be used in tau-leap
    if (!m_allowedMovements) {
        // if movement banned
        ++m_daysSinceMovementBan;
    }
}

bool Node::getAllowedToMove() const {
    return m_allowedMovements;
}

int Node::getDaysSinceMovementBan() const {
    return m_daysSinceMovementBan;
}

void Node::checkIsInfected() {
    /* Double check the infection status of the Node, and adjust as necessary.
     * Three scenarios:
     * - Not infected
     * - Infected with 1 serotype
     * - Infected with 2+ serotypes -> prefer one already infected with, otherwise choose randomly */

    std::vector<size_t> infectedIndices;
    for (size_t i = 0; i < m_serotypeCompartments.size(); ++i) {
        if (m_serotypeCompartments[i].getIsInfected()) {
            infectedIndices.push_back(i);
        }
    }
    auto numInfected{infectedIndices.size()};
    switch (numInfected) {
        case 0:
            m_isInfected = false;
            m_serotypeInfectedWith = "none";
            m_infectionDetectable = true;
            break;
        case 1:
            m_isInfected = true;
            m_serotypeInfectedWith = m_serotypeCompartments[infectedIndices[0]].getSerotype();
            break;
        default:
            // hopefully never gets here. Should lof if it does, once I implement logging.
            if (m_isInfected) {
                for (auto &scm : m_serotypeCompartments) {
                    if (scm.getSerotype() != m_serotypeInfectedWith) {
                        scm.uninfect();
                    }
                }
            } else {
                m_isInfected = true;
                std::uniform_int_distribution<size_t> dist(0, numInfected - 1);
                size_t selected{dist(*(m_pRNG))};
                for (size_t i = 0; i < numInfected; ++i) {
                    size_t ind{infectedIndices[i]};
                    if (i == selected) {
                        m_serotypeInfectedWith = m_serotypeCompartments[ind].getSerotype();
                    } else {
                        m_serotypeCompartments[ind].uninfect();
                    }
                }
            }
            break;
    }
}

bool Node::detectInfected() {
    m_hasBeenDetected = (m_isInfected && m_infectionDetectable &&
                         (m_daysSinceInfected > int(m_config->m_detectionDelay)));
    return m_hasBeenDetected;
}

void Node::uninfect() {
    if (m_isInfected) {
        auto it{findSerotype(m_serotypeInfectedWith)};
        (*it).uninfect();
    }
}

const std::string &Node::getSerotypeInfectedWith() const {
    return m_serotypeInfectedWith;
}

bool Node::canBeInfectedWithSerotype(const std::string &infectingSerotype) const {
    // ASSUMPTION: Already infected premises can still be infected with the same serotype, so long as S > 0
    // it just converts one S -> E
    // ASSUMPTIOM: Node infected with another serotype cannot be infected
    auto it{findSerotype(infectingSerotype)};
    bool susceptibleToSerotype = (*it).getCompartmentPopulation(Status::SUS) > 0;
    bool noConflictingSerotype = (m_serotypeInfectedWith == infectingSerotype || m_serotypeInfectedWith == "none");
    return ((susceptibleToSerotype && noConflictingSerotype));
}

Cell *Node::getCellIn() const {
    return m_parentCell;
}

void Node::setCellIn(Cell *cell) {
    m_parentCell = cell;
}

void Node::setGridIn(Grid *gridPtr) {
    m_parentGrid = gridPtr;
}

int Node::getNumCattle() const {
    return m_numCattle;
}

double Node::getCattleInfectionProb(const Node *otherNode) {
    // This is intended to calculate the probability of this node infecting the other node.
    if (m_isInfected) {
        // distance
        double distance = common::euclideanDistance(m_location, otherNode->getLocation());
        double kernel = m_parentGrid->powerLawKernel(distance);
        // transmission_pars take into account number of infectious cattle
        return common::infectionProbabilty(m_transmissionPars.first, m_transmissionPars.second, kernel);
    } else {
        return 0.0;
    }
}


void Node::setIsSentinel(bool newVal) {
    m_isSentinelNode = newVal;
}

void Node::updateDaysSinceVaccinated() {
    if (m_hasBeenVaccinated) {
        ++m_daysSinceVaccination;
    }
}

int Node::getDaysSinceVaccinated() const {
    return m_daysSinceVaccination;
}

bool Node::hasBeenVaccinated() const {
    return m_hasBeenVaccinated;
}

void Node::updateTransmission() {

    if (m_isInfected) {
        auto it{findSerotype(m_serotypeInfectedWith)};
        int totalInfectious = (*it).getSumInfectious();

        m_transmissionPars.first = totalInfectious *
                m_config->m_serotypeInfectionParameters.at(
                        m_serotypeInfectedWith).transmission;
        // don't multiply susceptibility, is multiplied by num susceptible in recipient node.
        m_transmissionPars.second = m_config->m_serotypeInfectionParameters.at(
                m_serotypeInfectedWith).susceptibility;
    } else {
        m_transmissionPars.first = 0.0;
        m_transmissionPars.second = 0.0;
    }

}

bool Node::isInfectiousWith(const std::string &serotype) const {
    auto it{findSerotype(serotype)};
    return (m_isInfected && (*it).getSumInfectious() > 0);
}

std::vector<int> Node::getSumCompartments(std::string_view serotype) const {
    auto it{findSerotype(serotype)};
    std::vector<int> compartmentSums{(*it).getCompartmentSums()};
    return compartmentSums;
}

bool Node::getIsQueuedForVaccination() const {
    return m_isQueuedForVaccination;
}

void Node::setIsQueuedForVaccination(bool mIsQueuedForVaccination) {
    m_isQueuedForVaccination = mIsQueuedForVaccination;
}

bool Node::getHasBeenDetected() const {
    return m_hasBeenDetected;
}

void Node::readState(std::vector<std::string> const &tokens) {
    // 25 headers (and position) are:
    // node_id (0), km_north (1), km_east (2), serotype (3), population (4), is_infected (5),
    // serotype_infected (6), cattle_infectiousness (7), cattle_susceptibility (8),
    // infection_detectable (9), infection_detected (10), days_infected (11), allowed_movement (12),
    // days_banned (13), vaccination_queue (14), been_vaccinated (15), days_vaccinated (16),
    // MAT (17), SUS (18), EXP (19), INF (20), REC (21), SUSVAC (22), RECVAC (23), CAR (24)
    m_isInfected = common::trueFalseToBool(tokens[5]);
    m_serotypeInfectedWith = tokens[6];
    m_transmissionPars.first = std::stod(tokens[7]);
    m_transmissionPars.second = std::stod(tokens[8]);
    m_infectionDetectable = common::trueFalseToBool(tokens[9]);
    m_hasBeenDetected = common::trueFalseToBool(tokens[10]);
    m_daysSinceInfected = std::stoi(tokens[11]);
    m_allowedMovements = common::trueFalseToBool(tokens[12]);
    m_daysSinceMovementBan = std::stoi(tokens[13]);
    m_isQueuedForVaccination = common::trueFalseToBool(tokens[14]);
    m_hasBeenVaccinated = common::trueFalseToBool(tokens[15]);
    m_daysSinceVaccination = std::stoi(tokens[16]);
}

YAML::Emitter &operator<<(YAML::Emitter &out, Node const &node) {
    out << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "id"; // node id
    out << YAML::Value << node.m_epiID;
    out << YAML::Key << "y_coord"; // y-coord
    out << YAML::Value << node.m_location.getYCoord();
    out << YAML::Key << "x_coord"; // x-coord
    out << YAML::Value << node.m_location.getXCoord();
    out << YAML::Key << "num_cattle"; // num cattle
    out << YAML::Value << node.m_numCattle;
    out << YAML::Key << "is_infected"; // is infected
    out << YAML::Value << node.m_isInfected;
    out << YAML::Key << "serotype_infected_with";
    out << YAML::Value << node.m_serotypeInfectedWith;
    out << YAML::Key << "transmission_first";
    out << YAML::Value << node.m_transmissionPars.first;
    out << YAML::Key << "transmission_second";
    out << YAML::Value << node.m_transmissionPars.second;
    out << YAML::Key << "infection_detectable";
    out << YAML::Value << node.m_infectionDetectable;
    out << YAML::Key << "has_been_detected";
    out << YAML::Value << node.m_hasBeenDetected;
    out << YAML::Key << "days_since_infected";
    out << YAML::Value << node.m_daysSinceInfected;
    out << YAML::Key << "allowed_movements";
    out << YAML::Value << node.m_allowedMovements;
    out << YAML::Key << "days_since_movement_ban";
    out << YAML::Value << node.m_daysSinceMovementBan;
    out << YAML::Key << "queued_for_vaccination";
    out << YAML::Value << node.m_isQueuedForVaccination;
    out << YAML::Key << "has_been_vaccinated";
    out << YAML::Value << node.m_hasBeenVaccinated;
    out << YAML::Key << "days_since_vaccination";
    out << YAML::Value << node.m_daysSinceVaccination;
    for (auto const &scm : node.m_serotypeCompartments) {
        out << YAML::Key << "serotype_" + scm.getSerotype();
        out << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "serotype";
        out << YAML::Value << scm.getSerotype();
        out << YAML::Key << "MAT";
        out << YAML::Value << scm.getSubCompartments(Status::MAT);
        out << YAML::Key << "SUS";
        out << YAML::Value << scm.getSubCompartments(Status::SUS);
        out << YAML::Key << "EXP";
        out << YAML::Value << scm.getSubCompartments(Status::EXP);
        out << YAML::Key << "INF";
        out << YAML::Value << scm.getSubCompartments(Status::INF);
        out << YAML::Key << "CAR";
        out << YAML::Value << scm.getSubCompartments(Status::CAR);
        out << YAML::Key << "REC";
        out << YAML::Value << scm.getSubCompartments(Status::REC);
        out << YAML::Key << "SUSVAC";
        out << YAML::Value << scm.getSubCompartments(Status::SUSVAC);
        out << YAML::Key << "RECVAC";
        out << YAML::Value << scm.getSubCompartments(Status::RECVAC);
        out << YAML::EndMap;
    }
    out << YAML::EndMap;
    return out;
}

