//
// Created by glen on 03/02/2020.
//

#include <algorithm>
#include "Common.h"
#include "Cell.h"
#include "Grid.h"
#include "Config.h"


Cell::Cell(Point lower_left, double width, double height, Grid *parentGrid, const Config *config, std::mt19937 *rng)
        : m_cellID(std::to_string(s_cellCount++)), m_config(config),
          m_parentGrid(parentGrid), m_pRNG(rng), m_cellWidth(width), m_cellHeight(height) {

    m_cellCorners.push_back(lower_left);
    m_cellCorners.emplace_back(lower_left.getXCoord(), lower_left.getYCoord() + height); // upper left
    m_cellCorners.emplace_back(lower_left.getXCoord() + width, lower_left.getYCoord() + height); // upper right
    m_cellCorners.emplace_back(lower_left.getXCoord() + width, lower_left.getYCoord()); // lower right
    m_boundsX.emplace_back(lower_left.getXCoord());
    m_boundsX.emplace_back(lower_left.getXCoord() + width);
    m_boundsY.emplace_back(lower_left.getYCoord());
    m_boundsY.emplace_back(lower_left.getYCoord() + height);
    for (const auto &s : m_config->m_serotypes) {
        m_serotypeSusceptibleNodes.emplace(s, std::vector<Node *>());
        m_serotypeInfectedNodes.emplace(s, std::vector<Node *>());
    }
}

void Cell::addNode(Node *newNode) {
    m_allNodes.push_back(newNode);
    newNode->setCellIn(this);
    ++m_numNodes;
    assignNode(newNode); // assign to susceptible/infectious nodes
}

void Cell::addNodes(std::vector<Node> &newNodes) {
    m_allNodes.reserve(m_allNodes.size() + newNodes.size());
    for (auto &node : newNodes) {
        m_allNodes.push_back(&node);
        node.setCellIn(this);
        assignNode(&node);
    }
    m_numNodes += newNodes.size();
}

void Cell::addNodes(std::vector<Node *> &newNodes) {
    m_allNodes.insert(m_allNodes.end(), newNodes.begin(), newNodes.end());
    for (auto &nodePtr : newNodes) {
        nodePtr->setCellIn(this);
        assignNode(nodePtr);
    }
    m_numNodes += newNodes.size();
}

void Cell::assignNode(Node *nodePtr) {
    // assign to susceptible/infectious, and also update max_susceptibility/infectiousness while at it
    if (nodePtr->isInfected()) {
        m_allInfectiousNodes.push_back(nodePtr);
        ++m_numInfectiousNodes;
        // update serotypes present in the Cell
        if (std::find(m_serotypesPresent.begin(), m_serotypesPresent.end(), nodePtr->getSerotypeInfectedWith()) ==
            m_serotypesPresent.end()) {
            m_serotypesPresent.push_back(nodePtr->getSerotypeInfectedWith());
        }
    } else {
        m_susceptibleNodes.push_back(nodePtr);
        ++m_numSusceptibleNodes;
    }
    // check if need to update maxPossibleNodeSize for calculation of max
    if (nodePtr->getNumCattle() > m_maxPossibleNodeSize) {
        if (m_sentinelNode != nullptr) {
            m_sentinelNode->setIsSentinel(false);
        }
        m_sentinelNode = nodePtr;
        nodePtr->setIsSentinel(true);
        m_maxPossibleNodeSize = nodePtr->getNumCattle();
        updateMaxParameters();
    }
}

void Cell::updateMaxParameters() {
    m_maxParameters[0] = m_config->m_maxTransmission * m_maxPossibleNodeSize;
    m_maxParameters[1] = m_config->m_maxSusceptibility * m_maxPossibleNodeSize;
}

std::vector<Cell> Cell::subdivide(size_t minNodesPerCell) {
    /* Recursively subdivide this Cell until approximately around calculated optimal number of Nodes are
     * within each Cell.
     * - Generate subdivision Cells and calculate how many Nodes they would contain
     * - If closer to optimal value than parent Cell without going under, actually add the Nodes and call subdivide.
     * - If variance from optimal not less, return self?
     * - If number of nodes <= optimal value, return self? */
    std::vector<Cell> completedCells;
    if (m_allNodes.size() > minNodesPerCell) {
        std::vector<Cell> childCells;
        double halfWidth = m_cellWidth / 2;
        double halfHeight = m_cellHeight / 2;
        // lower left cell
        childCells.emplace_back(m_cellCorners[0], halfWidth, halfHeight, m_parentGrid, m_config, m_pRNG);
        // upper left cell
        childCells.emplace_back(Point{m_cellCorners[0].getXCoord(), m_cellCorners[0].getYCoord() + halfHeight},
                                halfWidth, halfHeight, m_parentGrid, m_config, m_pRNG);
        // upper right cell
        childCells.emplace_back(
                Point{m_cellCorners[0].getXCoord() + halfWidth, m_cellCorners[0].getYCoord() + halfHeight},
                halfWidth, halfHeight, m_parentGrid, m_config, m_pRNG);
        // lower cell
        childCells.emplace_back(Point{m_cellCorners[0].getXCoord() + halfWidth, m_cellCorners[0].getYCoord()},
                                halfWidth, halfHeight, m_parentGrid, m_config, m_pRNG);

        // calculate variance etc.
        std::vector<std::vector<Node *>> cellNodes(4);
        std::vector<double> logCSizes; // log of cell size
        std::vector<double> cSizes; // num nodes in each cell
        for (size_t indexCell = 0; indexCell < 4; ++indexCell) {
            for (const auto &nodePtr : m_allNodes) {
                if (childCells[indexCell].nodeIsWithinCell(nodePtr)) {
                    cellNodes[indexCell].push_back(nodePtr);
                }
            }
            // count number of nodes that end up in each cell (if > 0)
            size_t numNodesInCell{cellNodes[indexCell].size()};
            if (numNodesInCell > 0) {
                logCSizes.push_back(std::log(numNodesInCell));
                cSizes.push_back(double(numNodesInCell));
            }
        }

        // calculate new + original variance from target value
        double cellVar = common::variance_from(logCSizes, std::log(minNodesPerCell));
        double parentVar = std::pow((std::log(m_numNodes) - std::log(minNodesPerCell)), 2);
        if (cellVar < parentVar) {
            // actually add nodes
            std::vector<Cell> childrenWithNodes;
            for (size_t indexCell = 0; indexCell < 4; ++indexCell) {
                if (not cellNodes[indexCell].empty()) {
                    childCells[indexCell].addNodes(cellNodes[indexCell]);
                    childrenWithNodes.push_back(childCells[indexCell]);
                }
            }

            // subdivide again
            for (auto &child : childrenWithNodes) {
                std::vector<Cell> grandchildren = child.subdivide(minNodesPerCell);
                completedCells.insert(completedCells.end(), grandchildren.begin(), grandchildren.end());
            }
        } else {
            completedCells.insert(completedCells.end(), *this);
        }
    } else {
        completedCells.insert(completedCells.end(), *this);
    }
    return completedCells;
}

bool Cell::nodeIsWithinCell(const Node *nodeToCheck) const {
    return (pointIsWithin(nodeToCheck->getLocation()));
}

bool Cell::pointIsWithin(const Point &pointToCheck) const {
    return (pointIsWithinVerticalRange(pointToCheck) && pointIsWithinHorizontalRange(pointToCheck));
}

bool Cell::pointIsWithinVerticalRange(const Point &pointToCheck) const {
    return (pointToCheck.getXCoord() >= m_boundsX[0] && pointToCheck.getXCoord() <= m_boundsX[1]);
}

bool Cell::pointIsWithinHorizontalRange(const Point &pointToCheck) const {
    return (pointToCheck.getYCoord() >= m_boundsY[0] && pointToCheck.getYCoord() <= m_boundsY[1]);
}

double Cell::calculateShortestDistance(const Cell &otherCell) {
    // calculates the shortest distance inside a grid of heterogenously sized cells
    bool shortestDistFound{false};
    std::vector<double> newDistances;
    newDistances.reserve(8);

    std::vector<Point> otherCorners{otherCell.getCorners()};
    for (const auto &otherPoint : otherCorners) {
        for (const auto &ownPoint : m_cellCorners) {
            // if share a point then adjacent, 0.0 distance
            if (otherPoint == ownPoint) {
                newDistances.push_back(0.0);
                shortestDistFound = true;
                break;
            }
        }
    }
    // check for horizontal alignment, in which case shortest distance is straight line between edges
    if (!shortestDistFound) {
        for (const auto &otherPoint : otherCorners) {
            if (pointIsWithinHorizontalRange(otherPoint)) {
                newDistances.push_back(std::abs(m_cellCorners[0].getYCoord() - otherPoint.getYCoord())); // lower left
                newDistances.push_back(std::abs(m_cellCorners[1].getYCoord() - otherPoint.getYCoord())); // upper left
                shortestDistFound = true;
            }
        }
        for (const auto &ownPoint : m_cellCorners) {
            if (otherCell.pointIsWithinHorizontalRange(ownPoint)) {
                newDistances.push_back(std::abs(otherCorners[0].getYCoord() - ownPoint.getYCoord())); // lower left
                newDistances.push_back(std::abs(otherCorners[1].getYCoord() - ownPoint.getYCoord())); // upper left
                shortestDistFound = true;
            }
        }
    }
    // check vertical alignment, although only if no horizontal alignment found as can't be both
    if (!shortestDistFound) {
        for (const auto &otherPoint : otherCorners) {
            if (pointIsWithinVerticalRange(otherPoint)) {
                newDistances.push_back(std::abs(m_cellCorners[0].getXCoord() - otherPoint.getXCoord())); // lower left
                newDistances.push_back(std::abs(m_cellCorners[2].getXCoord() - otherPoint.getXCoord())); // lower right
                shortestDistFound = true;
            }
        }
        for (const auto &ownPoint : m_cellCorners) {
            if (otherCell.pointIsWithinVerticalRange(ownPoint)) {
                newDistances.push_back(std::abs(otherCorners[0].getXCoord() - ownPoint.getXCoord())); // lower left
                newDistances.push_back(std::abs(otherCorners[2].getXCoord() - ownPoint.getXCoord())); // lower right
                shortestDistFound = true;
            }
        }
    }

    // no vertical or horizontal alignment, need to use euclidean distance
    if (!shortestDistFound) {
        for (const auto &ownPoint : m_cellCorners) {
            for (const auto &otherPoint : otherCorners) {
                newDistances.push_back(common::euclideanDistance(ownPoint, otherPoint));
            }
        }
    }
    return (*std::min_element(newDistances.begin(), newDistances.end()));
}

std::vector<Point> Cell::getCorners() const {
    return m_cellCorners;
}

double Cell::getMaxInfectiousness() const {
    return m_maxParameters[0];
}

double Cell::getMaxSusceptibility() const {
    return m_maxParameters[1];
}

size_t Cell::calculateDistanceIndex() const {
    for (size_t i = 0; i < m_parentGrid->m_allCells.size(); ++i) {
        if (m_parentGrid->m_allCells[i].getID() == m_cellID) {
            return i;
        }
    }
    // if hasn't found itself (which should never happen), cannot get shared_ptr to itself so throw exception
    throw std::runtime_error("Cell object called which is not owned by a Grid.");
}

std::string Cell::getID() const {
    return m_cellID;
}

size_t Cell::getDistanceIndex() const {
    if (m_hasCalculatedDistanceIndex) {
        return m_distanceIndex;
    } else {
        m_distanceIndex = calculateDistanceIndex();
        return m_distanceIndex;
    }
}

double Cell::getShortestDistanceToCell(const Cell &otherCell) const {
    // use this->getDistanceIndex() to make sure has value, as caches
    return m_parentGrid->m_matrixDistanceCellToCell(this->getDistanceIndex(), otherCell.getDistanceIndex());
}

void Cell::setGridIn(Grid *gridPtr) {
    m_parentGrid = gridPtr;
}

void Cell::reassignNodes(std::vector<Node> &newNodes) {
    // when copying a Cell as part of copying a Grid object, the new Cell needs to point to the Nodes actually
    // in that Grid object, which are not the same ones as those in the copied Grid object.
    // So reset, and reallocate Nodes.
    m_allNodes.clear();
    m_allInfectiousNodes.clear();
    m_susceptibleNodes.clear();
    m_numNodes = 0;
    m_numSusceptibleNodes = 0;
    m_numInfectiousNodes = 0;
    m_serotypesPresent.clear();
    for (auto &node : newNodes) {
        if (nodeIsWithinCell(&node)) {
            addNode(&node);
        }
    }
}

void Cell::updateNodes() {
    // reassign to susceptible/infectious vectors anew, and recalculate transmission parameters
    m_susceptibleNodes.clear();
    m_allInfectiousNodes.clear();
    m_numSusceptibleNodes = 0;
    m_numInfectiousNodes = 0;
    m_serotypesPresent.clear();
    for (const auto &node : m_allNodes) {
        assignNode(node);
    }
}

bool Cell::hasInfectious() const {
    return m_numInfectiousNodes > 0;
}

double Cell::getOverestimatedP(const Cell &otherCell) const {
    return m_parentGrid->m_matrixOverestimatedProbabilities(getDistanceIndex(), otherCell.getDistanceIndex());
}

std::vector<std::string> Cell::getSerotypesPresent() const {
    return m_serotypesPresent;
}

void Cell::addSerotypePresent(const std::string &serotype) {
    if (std::find(m_serotypesPresent.begin(), m_serotypesPresent.end(), serotype) == m_serotypesPresent.end()) {
        m_serotypesPresent.push_back(serotype);
    }
}

int Cell::getMaxNodeSize() const {
    return m_maxPossibleNodeSize;
}

void Cell::updateMaxNodeSize(int newSize) {
    m_maxPossibleNodeSize = newSize;
    updateMaxParameters();
}

void
Cell::calculateSerotypeInfectionProbabilities(std::unordered_map<Node *, std::map<std::string, double>> &serotypeProbs,
                                              const Cell &infCell, double conditionalP) const {
    /* Purpose: to calculate the serotype-specific probability of infection for each possible serotype-susceptible node
     * This is necessary as infection needs to take into account multiple different serotypes.
     * ConditonalP should be 1.0 for within-cell local spread, and input during conditional subsampling for between-cells */

    for (const auto &sero : infCell.getSerotypesPresent()) {


        std::vector<Node *> susToSero{m_serotypeSusceptibleNodes.at(sero)};
        std::vector<Node *> infWithSero{infCell.m_serotypeInfectedNodes.at(sero)};
        int numSus = static_cast<int>(susToSero.size());

        // short circuit, only run calculations if there are susceptible
        if (numSus > 0 && not infWithSero.empty()) {
            std::vector<Node *> sampleNodes;
            if (conditionalP < 1.0) {
                std::binomial_distribution<int> dist(numSus, conditionalP);
                int numInfected = dist(*(m_pRNG));
                if (numInfected > 0) {
                    std::sample(susToSero.begin(), susToSero.end(), std::back_inserter(sampleNodes), numInfected,
                                *(m_pRNG));
                }
            } else {
                // shortcut, when conditionalP is 1.0 will sample all anyway
                sampleNodes = std::move(susToSero);
            }
            for (const auto &sampleNode : sampleNodes) {
                /* 3 possibilities:
                 * 1 - node map is present, as is serotype specific probability
                 * 2 - node map is present, but hasn't assessed specific serotype yet
                 * 3 - node map is not present*/
                auto nodePMap{serotypeProbs.find(sampleNode)};
                if(nodePMap != serotypeProbs.end()) {
                    auto current{nodePMap->second.find(sero)};
                    if(current != nodePMap->second.end()) {
                        double negPreviousP{1.0 - current->second}; // convert back to negative prob
                        double currentP{1.0};
                        for (const auto &infNode : infWithSero) {
                            double infProb = infNode->getCattleInfectionProb(sampleNode);
                            currentP *= (1.0 - infProb);
                        }
                        currentP = (1.0 - currentP) / conditionalP; // positive prob
                        double negCurrentP = 1.0 - currentP;
                        double totalP{negCurrentP * negPreviousP};
                        current->second = 1.0 - totalP;
                    } else {
                        // this specific serotype hasn't been assessed yet
                        double initialP{1.0};
                        for (const auto &infNode : infWithSero) {
                            double infProb = infNode->getCattleInfectionProb(sampleNode);
                            // already made sure same serotype as sero
                            initialP *= (1.0 - infProb);
                        }
                        nodePMap->second[sero] = (1.0 - initialP) / conditionalP;
                    }
                } else {
                    std::map<std::string, double> newNodePMap;
                    double initialP{1.0};
                    // calculate probability at least 1 recipient cattle _isn't_ infected
                    for (const auto &infNode : infWithSero) {
                        double infProb = infNode->getCattleInfectionProb(sampleNode);
                        // already made sure same serotype as sero
                        initialP *= (1.0 - infProb);
                    }
                    newNodePMap[sero] = (1.0 - initialP) / conditionalP;
                    serotypeProbs[sampleNode] = newNodePMap;
                }
            }
        }
    }
}

void Cell::updateSerotypeNodeMap() {
    // sort nodes into serotype-specific susceptible/infected maps
    for (const auto &s : m_config->m_serotypes) {
        m_serotypeSusceptibleNodes.at(s).clear();
        std::copy_if(m_allNodes.begin(), m_allNodes.end(), std::back_inserter(m_serotypeSusceptibleNodes.at(s)),
                     [&](const Node *node) { return node->canBeInfectedWithSerotype(s); });
        // only update infectious nodes if there are ones of this serotype present
        if (std::find(m_serotypesPresent.begin(), m_serotypesPresent.end(), s) != m_serotypesPresent.end()) {
            m_serotypeInfectedNodes.at(s).clear();
            std::copy_if(m_allInfectiousNodes.begin(), m_allInfectiousNodes.end(),
                         std::back_inserter(m_serotypeInfectedNodes.at(s)),
                         [&](const Node *node) { return node->isInfectiousWith(s); });
        }
    }
}

bool Cell::operator==(const Cell &rhs) const {
    return m_cellID == rhs.m_cellID;
}

bool Cell::operator!=(const Cell &rhs) const {
    return !(rhs == *this);
}



