//
// Created by glen on 03/02/2020.
//

#include <iomanip>
#include <sstream>
#include <algorithm>
#include "Grid.h"
#include "Common.h"
#include "Config.h"

/// Constructor
Grid::Grid(const Config *config)
        : m_config(config), m_reporting(m_config->m_reporting.get()), m_pRNG(m_config->m_pRNG.get()) {
    m_currentDay = m_config->m_startDay;
    readNodeFile();
    calculateDimensions();
    calculateOptimalNodesPerCell();
    createDynamicCells();
    calculateShortestDistances();
    calculateKernelValues();
    initNodeRadii();
}

/// Copy Constructor
Grid::Grid(const Grid &rhs) : m_config(rhs.m_config), m_reporting(rhs.m_reporting), m_pRNG(rhs.m_pRNG) {
    // copy all the members
    m_currentDay = rhs.m_currentDay;
    m_verticalLength = rhs.m_verticalLength;
    m_horizontalLength = rhs.m_horizontalLength;
    m_longestSide = rhs.m_longestSide;
    m_longestDistance = rhs.m_longestDistance;
    m_numCells = rhs.m_numCells;
    m_optimalNodesPerCell = rhs.m_optimalNodesPerCell;
    m_boundsX = rhs.m_boundsX;
    m_boundsY = rhs.m_boundsY;
    // don't copy m_pointsX/Y as no use once m_boundsX/Y calculated
    m_lowerLeft = rhs.m_lowerLeft;
    m_upperLeft = rhs.m_upperLeft;
    m_upperRight = rhs.m_upperRight;
    m_lowerRight = rhs.m_lowerRight;

    m_matrixDistanceCellToCell = rhs.m_matrixDistanceCellToCell;
    m_matrixOverestimatedProbabilities = rhs.m_matrixOverestimatedProbabilities;

    m_allNodes.reserve(rhs.m_allNodes.size());
    for (const auto &node : rhs.m_allNodes) {
        m_allNodes.emplace_back(node); // uses default copy constructor
        m_allNodes.back().setGridIn(this);
        if (m_allNodes.back().isInfected()) {
            m_infectiousNodes.push_back(&m_allNodes.back());
        }
        if (node.getIsQueuedForVaccination()) {
            m_vaccinationQueue.push(&m_allNodes.back());
        }
        m_nodeIDLookup[node.getID()] = &m_allNodes.back();
    }

    // Cells own viewing pointers, so need to update to point new memory (Grid, Nodes).
    m_allCells.reserve(rhs.m_allCells.size());
    for (size_t i = 0; i < rhs.m_allCells.size(); ++i) {
        m_allCells.emplace_back(rhs.m_allCells[i]);
        m_allCells[i].setGridIn(this);
        m_allCells[i].calculateDistanceIndex();
        m_allCells[i].reassignNodes(m_allNodes);
    }
}

/// Copy Assignment
Grid &Grid::operator=(const Grid &rhs) {
    // Self-assignment guard
    if (this == &rhs) {
        return *this;
    }
    // copy all the members
    m_config = rhs.m_config;
    m_reporting = rhs.m_reporting;
    m_pRNG = rhs.m_pRNG;
    m_currentDay = rhs.m_currentDay;
    m_verticalLength = rhs.m_verticalLength;
    m_horizontalLength = rhs.m_horizontalLength;
    m_longestSide = rhs.m_longestSide;
    m_longestDistance = rhs.m_longestDistance;
    m_numCells = rhs.m_numCells;
    m_optimalNodesPerCell = rhs.m_optimalNodesPerCell;
    m_boundsX = rhs.m_boundsX;
    m_boundsY = rhs.m_boundsY;
    // don't copy m_pointsX/Y as no use once m_boundsX/Y calculated
    m_lowerLeft = rhs.m_lowerLeft;
    m_upperLeft = rhs.m_upperLeft;
    m_upperRight = rhs.m_upperRight;
    m_lowerRight = rhs.m_lowerRight;

    m_allNodes.reserve(rhs.m_allNodes.size());
    for (const auto &node : rhs.m_allNodes) {
        m_allNodes.emplace_back(node); // uses default copy constructor
        m_allNodes.back().setGridIn(this);
        if (m_allNodes.back().isInfected()) {
            m_infectiousNodes.push_back(&m_allNodes.back());
        }
        if (node.getIsQueuedForVaccination()) {
            m_vaccinationQueue.push(&m_allNodes.back());
        }
        m_nodeIDLookup[node.getID()] = &m_allNodes.back();
    }

    m_matrixDistanceCellToCell = rhs.m_matrixDistanceCellToCell;
    m_matrixOverestimatedProbabilities = rhs.m_matrixOverestimatedProbabilities;

    // // Cells own viewing pointers, so need to update to point new memory (Grid, Nodes).
    m_allCells.reserve(rhs.m_allCells.size());
    for (size_t i = 0; i < rhs.m_allCells.size(); ++i) {
        m_allCells.emplace_back(rhs.m_allCells[i]);
        m_allCells[i].setGridIn(this);
        m_allCells[i].calculateDistanceIndex();
        m_allCells[i].reassignNodes(m_allNodes);
    }

    return *this;
}

void Grid::readNodeFile() {
    if (m_config->m_savedNodeInput) {
        readSavedNodeFile();
    } else {
        readCleanNodeFile();
    }
}

void Grid::readCleanNodeFile() {
    if (m_config->m_verbosity >= 1) {
        std::cout << "Reading node file..." << std::endl;
    }
    std::shared_ptr<Grid> tempPtr = std::make_shared<Grid>((*this));
    std::ifstream file(m_config->m_nodeFile);
    if (file) {
        std::string line;
        while (getline(file, line)) {
            if (!std::isalpha(line[0])) {
                std::stringstream rowstream;
                rowstream << line;
                std::vector<std::string> tokens{common::split(rowstream, ',')};
                std::string ID{common::trimEnds(tokens[0])};
                double y{std::stod(tokens[1])};
                double x{std::stod(tokens[2])};
                m_pointsX.push_back(x);
                m_pointsY.push_back(y);
                Point location{x, y};
                int numCattle = std::stoi(tokens[3]);
                m_allNodes.emplace_back(ID, location, numCattle, this, m_config, m_pRNG);
            }
        }
    } else {
        throw std::runtime_error(" Cannot open file: " + m_config->m_nodeFile);
    }
    // inserting into vector reallocates, so need to set map when memory fully set.
    for (auto &Node : m_allNodes) {
        m_nodeIDLookup[Node.getID()] = &Node;
    }
}

void Grid::readSavedNodeFile() {
    if (m_config->m_verbosity >= 1) {
        std::cout << "Reading saved node file..." << std::endl;
    }
    std::shared_ptr<Grid> tempPtr = std::make_shared<Grid>((*this));
    YAML::Node yamlNodes = YAML::LoadFile(m_config->m_nodeFile);
    int count{1};
    while (yamlNodes["node_" + std::to_string(count)]) {
        YAML::Node current;
        current.reset(yamlNodes["node_" + std::to_string(count)]);
        m_allNodes.emplace_back(current, this, m_config, m_pRNG);
        m_pointsX.push_back(current["x_coord"].as<double>());
        m_pointsY.push_back(current["y_coord"].as<double>());
        ++count;
    }
    // inserting into vector reallocates, so need to set map when memory fully set.
    for (auto &Node : m_allNodes) {
        m_nodeIDLookup[Node.getID()] = &Node;
    }

    // update tracking vectors/queues
    detectedNodes();
    infectedNodes();
    vaccinationQueueNodes();
}

void Grid::calculateOptimalNodesPerCell() {
    if (m_config->m_verbosity >= 1) {
        std::cout << "Calculating optimal number of nodes per cell..." << std::endl;
    }
    // Purpose: Given the number of nodes read in, what is the number of Cells that minimises the kernel calls?
    // ASSUMPTION/HEURISTIC: There is only one infected Node, and every other node is susceptible
    auto totalNumNodes{static_cast<double>(m_allNodes.size())};
    auto maxOptimal = static_cast<size_t>(std::sqrt(
            totalNumNodes)); // looking at one side here, will want >1 Node per Cell
    double minimumKernelCallsSoFar =
            static_cast<double>(totalNumNodes) - 1.0; // expected kernel calls when cells = 1, so all pairwise
    double optimalSoFar{totalNumNodes};
    int numIncreased{
            0}; // tracker, when on the upswing of optimal nodes, stop running as will only continue to increase from there
    for (size_t CellsOneSide = 2; CellsOneSide <= maxOptimal && numIncreased < 3; ++CellsOneSide) {
        auto nodesPerCell = totalNumNodes / (static_cast<double>(CellsOneSide) * static_cast<double>(CellsOneSide));
        double averageKernelCalls{calculateKernelCalls(CellsOneSide)};
        if (averageKernelCalls < minimumKernelCallsSoFar) {
            minimumKernelCallsSoFar = averageKernelCalls;
            optimalSoFar = nodesPerCell;
        } else if (averageKernelCalls > minimumKernelCallsSoFar) {
            numIncreased += 1;
        }
    }
    m_optimalNodesPerCell = static_cast<size_t>(optimalSoFar);
    if (m_config->m_verbosity >= 2) {
        std::cout << "Calculated Optimal Nodes / Cell: " << m_optimalNodesPerCell << '\n';
    }
}

double Grid::calculateKernelCalls(size_t CellsOnOneSide) {
    auto totalNodes{static_cast<double>(m_allNodes.size())};
    if (CellsOnOneSide == 1) {
        return totalNodes - 1.0;
    } else { // more than 1 cell
        std::vector<double> kernelCallsPerA;
        double nodesPerCell{totalNodes / (static_cast<double>(CellsOnOneSide) * static_cast<double>(CellsOnOneSide))};
        std::vector<size_t> coordsCells(CellsOnOneSide); // symmetrical grid, so can use coordinates for x and y
        std::iota(coordsCells.begin(), coordsCells.end(), 0);
        // cell A coordinates
        for (const auto &xA : coordsCells) {
            for (const auto &yA : coordsCells) {
                double sumKernelCallsA{0.0};
                sumKernelCallsA += nodesPerCell - 1; // pairwise within cell
                // cell B coordinates
                for (const auto &xB : coordsCells) {
                    for (const auto &yB : coordsCells) {
                        // calculate expected number of kernel calls between cell(xA, yA) and cell(xB, yB)
                        if (!(xA == xB and yA == yB)) { // not if cell is itself
                            double distanceBetweenCells{optimisationDistanceFunction(xA, yA, xB, yB, m_longestSide,
                                                                                     CellsOnOneSide)};
                            double kernelValue{powerLawKernel(distanceBetweenCells)};
                            double exponent =
                                    -m_config->m_maxTransmission * m_config->m_maxSusceptibility * kernelValue;
                            // use calculated max transmission/susceptibility for overestimate
                            double overestimatedP{common::oneMinusExp(exponent)};
                            sumKernelCallsA += (nodesPerCell * overestimatedP +
                                                1); // estimated number of kernel calls for this pair of cells
                        }
                    }
                }
                // record value for this cell A
                kernelCallsPerA.push_back(sumKernelCallsA);
            }
        }
        // calculate average kernel calls per A, and return
        return common::average(kernelCallsPerA);
    }
}

double Grid::powerLawKernel(double distance) {
    return (m_config->m_a / (1 + pow((distance / m_config->m_scale), m_config->m_shape)));
}

double Grid::optimisationDistanceFunction(const size_t &xA, const size_t &yA, const size_t &xB, const size_t &yB,
                                          double longestSide, size_t cellsOneSide) {
    size_t a = 0;
    size_t b = 0;

    if (xA == xB)
        a = 0;
    if (xA < xB)
        a = 1;
    else if (xA > xB)
        a = static_cast<size_t>(-1);

    if (yA == yB)
        b = 0;
    else if (yA < yB)
        b = 1;
    else if (yA > yB)
        b = static_cast<size_t>(-1);

    double t1 = ((static_cast<double>(xA + a)) * longestSide - (double) xB * longestSide) / double(cellsOneSide);
    double t2 = ((static_cast<double>(yA + b)) * longestSide - (double) yB * longestSide) / double(cellsOneSide);
    double d_sq = t1 * t1 + t2 * t2;
    double d = std::sqrt(d_sq);

    return d;

}

void Grid::calculateDimensions() {
    double minX{(*std::min_element(m_pointsX.begin(), m_pointsX.end()))};
    double minY{(*std::min_element(m_pointsY.begin(), m_pointsY.end()))};
    double maxX{(*std::max_element(m_pointsX.begin(), m_pointsX.end()))};
    double maxY{(*std::max_element(m_pointsY.begin(), m_pointsY.end()))};
    m_boundsX.push_back(minX * 0.995); // add/subtract so Nodes not at boundaries
    m_boundsX.push_back(maxX * 1.005); // minimum, then maximum
    m_boundsY.push_back(minY * 0.995);
    m_boundsY.push_back(maxY * 1.005);
    m_horizontalLength = m_boundsX[1] - m_boundsX[0]; // max - min
    m_verticalLength = m_boundsY[1] - m_boundsY[0];
    m_longestSide = std::max(m_horizontalLength, m_verticalLength);
    m_longestDistance = std::sqrt(std::pow(m_horizontalLength, 2) + std::pow(m_verticalLength, 2)); // diagonal
    m_lowerLeft = Point{m_boundsX[0], m_boundsY[0]};
    m_upperLeft = Point{m_boundsX[0], m_boundsY[1]};
    m_upperRight = Point{m_boundsX[1], m_boundsY[1]};
    m_lowerRight = Point{m_boundsX[1], m_boundsY[0]};

    m_pointsX.clear();
    m_pointsY.clear(); // no use for possibly thousands of doubles now
}

void Grid::createDynamicCells() {
    if (m_config->m_verbosity >= 1) {
        std::cout << "Creating cell layout..." << std::endl;
    }
    Cell supercell{m_lowerLeft, m_horizontalLength, m_verticalLength, this, m_config, m_pRNG};
    supercell.addNodes(m_allNodes);
    m_allCells = supercell.subdivide(m_optimalNodesPerCell);
    // vector has been reallocated, so memory is different
    for (auto &cell : m_allCells) {
        cell.reassignNodes(m_allNodes);
    }
    m_numCells = m_allCells.size();
    if (m_config->m_verbosity >= 2) {
        std::cout << m_numCells << " Cells created.\n";
    }
}

void Grid::calculateShortestDistances() {
    m_matrixDistanceCellToCell = Matrix<unsigned int>(0, m_numCells, m_numCells);
    for (size_t indexA = 0; indexA < m_numCells; ++indexA) {
        for (size_t indexB = 0; indexB < m_numCells; ++indexB) {
            auto cellA(m_allCells[indexA]);
            auto cellB(m_allCells[indexB]);
            double distance = cellA.calculateShortestDistance(cellB);
            m_matrixDistanceCellToCell(indexA, indexB) = static_cast<unsigned int>(distance);
            m_matrixDistanceCellToCell(indexB, indexA) = static_cast<unsigned int>(distance);
        }
    }
}

void Grid::calculateKernelValues() {
    m_matrixOverestimatedProbabilities = Matrix<double>(0.0, m_numCells, m_numCells);

    for (size_t indexInf = 0; indexInf < m_numCells; ++indexInf) {
        for (size_t indexSus = 0; indexSus < m_numCells; ++indexSus) {
            auto cellInf(m_allCells[indexInf]);
            auto cellSus(m_allCells[indexSus]);
            double overP = calculateOverestimatedP(cellInf, cellSus);
            m_matrixOverestimatedProbabilities(indexInf, indexSus) = overP;
            m_matrixOverestimatedProbabilities(indexSus, indexInf) = overP;
        }
    }
}

double Grid::calculateOverestimatedP(const Cell &cellInf, const Cell &cellSus) {
    // cell to cell version
    double distance = cellInf.getShortestDistanceToCell(cellSus);
    double kernelVal{powerLawKernel(distance)};
    double inf_j = cellInf.getMaxInfectiousness();
    double sus_i = cellSus.getMaxSusceptibility();
    double exponent = inf_j * sus_i * kernelVal * -1;

    return common::oneMinusExp(exponent);
}

void Grid::run(const std::string &name, bool isBurnIn) {

    if (m_config->m_verbosity >= 1) {
        std::time_t t = std::time(nullptr);
        std::cout << "\n\n[" << std::put_time(std::localtime(&t), "%H:%M:%S") << "] - Beginning: " << name << std::endl;
    }
    // adjust for whether burn-in is happening/has happened
    int endDay = m_config->m_endDay;
    if (isBurnIn) {
        endDay = m_config->m_startDay + m_config->m_burnIn.duration - 1;
    }

    // DAYS
    for (int day = m_currentDay; day <= endDay; ++day) {
        incrementDay(day, name);
        updateCells();
    }
    if (not isBurnIn || (isBurnIn && not m_config->m_burnIn.discard)) {
        m_currentDay = endDay;
    }
}

void Grid::incrementDay(int day, const std::string &runName) {
    std::string dayString{std::to_string(day)};

    if (m_config->m_verbosity >= 2) {
        std::cout << day << ", ";
    }

    // SEEDING
    for (const auto &seed : m_config->m_seedEvents) {
        if (seed.day == day) {
            seedInfection(seed, runName);
        }
    }

    // TAU LEAP - for each farm
    for (auto &node : m_allNodes) {
        node.tauLeap();
    }
    // LOCAL SPREAD
    runLocalSpread(runName, dayString);

    // MOVEMENTS
    if (m_config->m_runShipments) {
        runShipments(runName, day);
    }

    // FORCING
    if (m_config->m_forceInfection) {
        std::uniform_real_distribution<double> forcingDist(0.0, 1.0);
        double actual = forcingDist(*(m_pRNG));
        if (actual <= m_config->m_forcingRate) {
            infectionForcing(runName, dayString);
        }
    }

    // DETECTION
    detectedNodes();

    std::pair<int, int> vaccinationStats;
    std::pair<int, int> movementBanStats;
    if (day >= m_config->m_dayEffectiveFrom) {
        if (day == m_config->m_dayEffectiveFrom && m_config->m_verbosity >= 1) {
            std::cout << "\n--- Control Policies Now Implemented ---\n";
        }
        // REACTIVE CONTROL
        if (m_config->m_movementBansImplemented) {
            movementBanStats = runMovementBans();
        }
        if (m_config->m_reactiveVaccinationImplemented) {
            identifyReactiveVaccinationTargets();
        }

        // MASS VACCINATION
        if (m_config->m_massVaccinationImplemented &&
            ((static_cast<unsigned int>(day - m_config->m_startDay)) % m_config->m_massVaccinationInterval == 0) &&
            day != m_config->m_startDay) {
            identifyMassVaccinationTargets();
        }

        // VACCINATE TARGETS
        if (not m_vaccinationQueue.empty()) {
            vaccinationStats = runVaccination();
        }
    }

    // UPDATE
    // NB: daysSinceInfected increments up, so if just infected will now be 1
    incrementNodes();
    updateMovementBans(); // if > duration, unban
    infectedNodes(); // updates, since no more use of m_infectiousNodes this day

    // REPORTING - ACTUALISATION
    m_reporting->updateModelReport(m_allNodes, runName, dayString);
    m_reporting->addVaccinationSum(runName, dayString, vaccinationStats);
    m_reporting->addMovementBanStats(m_allNodes, runName, dayString, movementBanStats);
    if (m_config->m_extraOutputs.globalCompartmentSums || m_config->m_extraOutputs.nodeCompartmentSums) {
        m_reporting->updateCompartmentSums(m_allNodes, runName, dayString);
    }

    if (m_config->m_extraOutputs.detailedReports) {
        m_reporting->updateDetailedReport(runName, dayString, m_allNodes);
    }
}

void Grid::seedInfection(const SeedEvent &seed, const std::string &runName) {
    // TODO: more intelligent handling of serotype seeding, e.g. seed with multiple serotypes?

    int numEpiSeeded;

    if (seed.isFixed) {
        numEpiSeeded = seedFixedSeeds(seed, runName);
    } else {
        numEpiSeeded = seedRandomSeeds(seed, runName);
    }

    if (numEpiSeeded < 1 && m_infectiousNodes.empty()) {
        throw std::runtime_error("No nodes have been successfully seeded, model ending.");
    }
    std::endl(std::cout);
}

int Grid::seedFixedSeeds(const SeedEvent &seed, std::string const &runName) const {
    int numSeeded{0};
    if (m_config->m_verbosity >= 1) {
        std::cout << "\n--- Attempting to seed fixed node " << seed.fixedSeedID << ": ";
    }
    try {
        auto nodePtr = m_nodeIDLookup.at(seed.fixedSeedID);
        std::map<std::string, double> prob{{seed.serotype, 1.0}};
        auto result = nodePtr->attemptToInfect(prob);
        if (result.first) {
            nodePtr->getCellIn()->addSerotypePresent(seed.serotype);
            if (m_config->m_extraOutputs.infectionEvents) {
                m_reporting->addInfectionEvent(runName, std::to_string(seed.day), "fixed-seed",
                                               seed.serotype, "seed", nodePtr->getID(), false);
            }
            ++numSeeded;
            if (m_config->m_verbosity >= 2) {
                std::cout << "\n-----> Seeded node " << seed.fixedSeedID;
            }
        } else {
            std::cout << "\n-----> Could not seed node " << seed.fixedSeedID;
        }

    } catch (std::out_of_range &e) {
        std::cerr << "\nCould not find Node: " << seed.fixedSeedID
                  << " in epiunits. Node was not seeded with infection.";
    }
    return numSeeded;
}

void Grid::incrementNodes() {
    for (auto &node : m_allNodes) {
        if (node.getDaysSinceInfected() > static_cast<int>(m_config->m_maxNodeInfectionDuration)) {
            node.uninfect();
        }
        node.update();
    }
}

std::vector<ValidatedShipment> Grid::validateShipments(int day) {
    // find movements on same day, with units allowed to move
    std::vector<ValidatedShipment> validShipments;
    for (auto &record : *(m_config->m_shipmentRecords)) {
        if (record.day == day) { // check move exists
            auto sourcePtr{m_nodeIDLookup.find(record.sourceNodeID)}; // value_type is pair<string, shared_ptr<Node>>
            auto targetPtr{m_nodeIDLookup.find(record.targetNodeID)};
            bool episExist{sourcePtr != m_nodeIDLookup.end() && targetPtr != m_nodeIDLookup.end()};
            if (episExist && sourcePtr->second->getAllowedToMove() && targetPtr->second->getAllowedToMove()) {
                validShipments.emplace_back(&record, sourcePtr->second, targetPtr->second);
            }
        }
    }
    return validShipments;
}

void Grid::runShipments(const std::string &runName, int day) {
    auto validShipments{validateShipments(day)};
    for (auto &valid  : validShipments) {
        // validateMovements means no exceptions thrown or movement bans
        Shipment sourceShipment = valid.sourcePtr->generateShipment(valid.record->numMoved);
        ShipmentResult result = valid.targetPtr->acceptShipment(sourceShipment);

        // record if option selected
        if (result.targetNowInfected && m_config->m_extraOutputs.infectionEvents) {
            std::string infectionType{};
            if (sourceShipment.shipmentIsInfected && not result.fomiteTransmissionOccured) {
                infectionType = "shipment-animals";
            } else if (sourceShipment.shipmentIsInfected && result.fomiteTransmissionOccured) {
                infectionType = "shipment-both";
            } else if (not sourceShipment.shipmentIsInfected && result.fomiteTransmissionOccured) {
                infectionType = "shipment-fomite";
            }
            m_reporting->addInfectionEvent(runName, std::to_string(day), infectionType,
                                           sourceShipment.infectionSerotype,
                                           valid.record->sourceNodeID,
                                           valid.record->targetNodeID,
                                           result.targetInfectedPrior);
        }
    }
}

void Grid::detectedNodes() {
    std::vector<Node *> detected;
    for (auto &node : m_allNodes) {
        if (node.detectInfected()) {
            detected.push_back(&node);
        }
    }
    m_detectedNodes = detected;
}

void Grid::infectedNodes() {
    std::vector<Node *> infected;
    for (auto &node : m_allNodes) {
        if (node.isInfected()) {
            infected.push_back(&node);
        }
    }
    m_infectiousNodes = infected;
}

void Grid::identifyReactiveVaccinationTargets() {
    std::uniform_int_distribution<int> coverageDist(0, 100);
    for (auto &detectedNode : m_detectedNodes) {
        for (const auto &otherNodeID : Grid::s_vaccinationRadiusIDs.at(detectedNode->getID())) {
            auto nodePtr{m_nodeIDLookup.at(otherNodeID)};
            // Vaccinate if hasn't been detected or queued for vaccination, if has been vaccinated then time should be > earliest allowed time.
            bool vaccinate{not(
                    nodePtr->getHasBeenDetected() ||
                    nodePtr->getIsQueuedForVaccination() ||
                    (nodePtr->hasBeenVaccinated() &&
                     nodePtr->getDaysSinceVaccinated() < m_config->m_reactiveVaccinationMinInterval))
                           && coverageDist(*(m_pRNG)) < m_config->m_reactiveVaccinationCoverage
            };
            if (vaccinate) {
                nodePtr->setIsQueuedForVaccination(true);
                m_vaccinationQueue.push(nodePtr);
            }
        }
    }
}

std::pair<int, int> Grid::runVaccination() {
    int nodesVaccinated{0};
    int animalsVaccinated{0};
    for (int i = 0; i < m_config->m_dailyVaccinationCapacity && not m_vaccinationQueue.empty(); ++i) {
        ++nodesVaccinated;
        animalsVaccinated += m_vaccinationQueue.front()->vaccinate();
        m_vaccinationQueue.pop();
    }
    return std::pair(nodesVaccinated, animalsVaccinated);
}

std::pair<int, int> Grid::runMovementBans() {
    // purpose: (un)ban movement from known infected units according to configVars
    int totalBansOrdered{0};
    int totalBansComplied{0};
    std::uniform_int_distribution<int> complianceDist(0, 100);
    for (auto &detectedNode : m_detectedNodes) {
        detectedNode->setMovement(false);
        ++totalBansOrdered;
        ++totalBansComplied;
        for (const auto &nodeID : Grid::s_banRadiusIDs.at(detectedNode->getID())) {

            auto nodePtr{m_nodeIDLookup.at(nodeID)};
            // if not banned, or is banned but not today, ban again
            if (nodePtr->getAllowedToMove() ||
                (not nodePtr->getAllowedToMove() && nodePtr->getDaysSinceMovementBan() > 0)) {
                ++totalBansOrdered;
                bool actuallyComplies = (complianceDist(*(m_pRNG)) <= m_config->m_movementBansCompliance);
                if (actuallyComplies) {
                    ++totalBansComplied;
                    nodePtr->setMovement(false);
                }
            }
        }
    }
    return std::pair(totalBansOrdered, totalBansComplied);
}

void Grid::identifyMassVaccinationTargets() {
    std::uniform_int_distribution<int> dist(0, 100);
    for (auto &node : m_allNodes) {
        if ((not node.getIsQueuedForVaccination())) {
            int coverageDraw{dist(*(m_pRNG))};
            bool isActuallyVaccinated = (coverageDraw < m_config->m_massVaccinationCoverage);
            if (isActuallyVaccinated) {
                node.setIsQueuedForVaccination(true);
                m_vaccinationQueue.push(&node);
            }
        }
    }
}

void Grid::infectionForcing(const std::string &runName, const std::string &dayString) {
    // Purpose: this function infects a random susceptible Node with a random serotype.
    // main decides the logic of when it is called.
    std::uniform_int_distribution<size_t> dist(0, m_config->m_serotypes.size() - 1);
    std::string selectedSerotype{m_config->m_serotypes.at(dist(*(m_pRNG)))};
    auto iter = std::find_if(m_allNodes.begin(), m_allNodes.end(),
                             [&](const Node &node) {
                                 return (!node.isInfected());
                             });
    if (iter != m_allNodes.end()) {
        std::map<std::string, double> prob{{selectedSerotype, 1.0}};
        auto result = (*iter).attemptToInfect(prob);
        if (result.first && m_config->m_extraOutputs.infectionEvents) {
            m_reporting->addInfectionEvent(runName, dayString, "forcing",
                                           selectedSerotype, "0", (*iter).getID(),
                                           false);
        }
    }
}

void Grid::updateMovementBans() {
    for (auto &node : m_allNodes) {
        if (node.getDaysSinceMovementBan() > static_cast<int>(m_config->m_movementBanDuration)) {
            node.setMovement(true); // allowed movements
        }
    }
}

void Grid::updateCells() {
    for (auto &cell : m_allCells) {
        cell.updateNodes();
    }
}

void Grid::runLocalSpread(const std::string &runName, const std::string &dayName) {
    // uses conditional subsample algorithm from Sellman, 2018 ("Need for Speed")
    std::unordered_map<Node*, std::map<std::string, double>> nodeSerotypePs;
    // make sure to update serotype-specific vectors
    for (auto &cell : m_allCells) {
        cell.updateSerotypeNodeMap();
    }

    for (auto &infCell : m_allCells) {
        if (infCell.hasInfectious()) {
            for (auto &recipientCell : m_allCells) {
                if (!(infCell == recipientCell)) {
                    double overP = infCell.getOverestimatedP(recipientCell);
                    recipientCell.calculateSerotypeInfectionProbabilities(nodeSerotypePs,
                                                                          infCell, overP);
                } else {
                    recipientCell.calculateSerotypeInfectionProbabilities(nodeSerotypePs,
                                                                          recipientCell, 1.0);
                }
            }
        }
    }

    // after calculating all probabilities, run through and see which are infected
    for (auto&[nodePtr, nodePMap] : nodeSerotypePs) {
        bool infectedPrior = nodePtr->isInfected();
        auto result = nodePtr->attemptToInfect(nodePMap);
        if (result.first && m_config->m_extraOutputs.infectionEvents) {
            m_reporting->addInfectionEvent(runName, dayName, "local-spread-cs", result.second, "CS",
                                           nodePtr->getID(), infectedPrior);
        }
    }
}

void Grid::initNodeRadii() {
    for (auto &node1 : m_allNodes) {
        for (auto &node2 : m_allNodes) {
            double distance{common::euclideanDistance(node1.getLocation(), node2.getLocation())};
            if (distance <= m_config->m_reactiveVaccinationRadius) {
                Grid::s_vaccinationRadiusIDs[node1.getID()].emplace_back(node2.getID());
            }
            if (distance <= m_config->m_movementBansRadius) {
                Grid::s_banRadiusIDs[node1.getID()].emplace_back(node2.getID());
            }
        }
    }
}

int Grid::seedRandomSeeds(const SeedEvent &seed, std::string const &runName) {
    if (seed.numberRandom > static_cast<int>(m_allNodes.size())) {
        throw std::invalid_argument("number_seeds cannot exceed the number of nodes supplied in node_file.");
    }

    int numSeeded{0};
    if (m_config->m_verbosity >= 1) {
        std::cout << "\n--- Seeding random node(s): ";
    }

    std::map<std::string, double> prob{{seed.serotype, 1.0}};
    std::uniform_int_distribution<size_t> nodeIndexDist(0, m_allNodes.size());
    while (numSeeded < seed.numberRandom) {
        size_t randomIndex{nodeIndexDist(*(m_pRNG))};
        auto nodePtr = &m_allNodes.at(randomIndex);
        if (nodePtr->canBeInfectedWithSerotype(seed.serotype)) {
            auto result = nodePtr->attemptToInfect(prob);
            if (result.first) {
                nodePtr->getCellIn()->addSerotypePresent(seed.serotype);
                if (m_config->m_extraOutputs.infectionEvents) {
                    m_reporting->addInfectionEvent(runName, std::to_string(seed.day), "random-seed",
                                                   seed.serotype, "seed", nodePtr->getID(), false);
                }
                ++numSeeded;
                if (m_config->m_verbosity >= 1) {
                    std::cout << "\n-----> Seeded node " << nodePtr->getID();
                }
            } else {
                if (m_config->m_verbosity >= 1) {
                    std::cout << "\n-----> Could not seed node " << nodePtr->getID();
                }
            }
        }
    }
    return numSeeded;
}

void Grid::saveNodeState() {
    if (m_config->m_verbosity >= 1) {
        std::cout << "\n- Writing: Node Save State";
    }
    m_reporting->saveNodeState(m_allNodes);
}

void Grid::vaccinationQueueNodes() {
    std::queue<Node *> queued;
    for (auto &node : m_allNodes) {
        if (node.getIsQueuedForVaccination()) {
            queued.push(&node);
        }
    }
    m_vaccinationQueue = queued;
}








