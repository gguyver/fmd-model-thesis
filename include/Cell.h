//
// Created by glen on 03/02/2020.
//

#ifndef MODEL_SINGLESEROTYPE_CELL_H
#define MODEL_SINGLESEROTYPE_CELL_H

// FORWARD DECLARATIONS
class Grid;

class Config;

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include "Node.h"

class Cell {
private:
    static inline int s_cellCount = 0; // for cellID
    std::string m_cellID;
    const Config *m_config;
    Grid *m_parentGrid; // ptr to Grid the Cell is in
    std::mt19937 *m_pRNG;

    bool m_hasCalculatedDistanceIndex{false}; // to know whether need to calculate it or not
    mutable size_t m_distanceIndex{0}; // used to track where in Grid::m_allCells this Cell is

    // LOCATION/SIZE
    std::vector<Point> m_cellCorners; // lower_left, upper_left, upper_right, lower_right
    std::vector<double> m_boundsX, m_boundsY; // minimum then maximum coordinates in x and y planes
    double m_cellWidth{0.0};
    double m_cellHeight{0.0};


    // NODES
    std::vector<std::string> m_serotypesPresent;
    std::vector<Node *> m_allNodes;
    std::vector<Node *> m_susceptibleNodes;
    std::vector<Node *> m_allInfectiousNodes;
    // store nodes suceptible/infected with specific serotypes
    std::map<std::string, std::vector<Node *>> m_serotypeSusceptibleNodes;
    std::map<std::string, std::vector<Node *>> m_serotypeInfectedNodes;
    Node *m_sentinelNode{nullptr}; // the node which is taken as upper limit for number of cattle for overP
    size_t m_numNodes{0};
    size_t m_numSusceptibleNodes{0};
    size_t m_numInfectiousNodes{0};
    int m_maxPossibleNodeSize{0};
    std::vector<double> m_maxParameters{0.0,
                                        0.0}; // max possible infectiousness, then max possible susceptiblity, calculated using most infectious serotype

    void assignNode(Node *nodePtr);

    void updateMaxParameters();

public:

    /// Constructor
    Cell(Point lower_left, double width, double height, Grid *parentGrid, const Config *config, std::mt19937 *rng);

    void addNode(Node *newNode);

    void addNodes(std::vector<Node> &newNodes);

    void addNodes(std::vector<Node *> &newNodes);


    std::vector<Cell> subdivide(size_t minNodesPerCell); // recursively subdivides cells until minNodesPerCell reached.

    bool nodeIsWithinCell(const Node *nodeToCheck) const;

    bool pointIsWithin(const Point &pointToCheck) const;

    bool pointIsWithinVerticalRange(const Point &pointToCheck) const;

    bool pointIsWithinHorizontalRange(const Point &pointToCheck) const;

    double calculateShortestDistance(const Cell &otherCell);

    std::vector<Point> getCorners() const;

    int getMaxNodeSize() const;

    void updateMaxNodeSize(int);

    double getMaxInfectiousness() const;

    double getMaxSusceptibility() const;

    double getShortestDistanceToCell(const Cell &otherCell) const;

    bool hasInfectious() const;

    size_t calculateDistanceIndex() const;

    std::string getID() const;

    size_t getDistanceIndex() const;

    void setGridIn(Grid *gridPtr);

    void
    calculateSerotypeInfectionProbabilities(std::unordered_map<Node *, std::map<std::string, double>> &serotypeProbs,
                                            const Cell &infCell, double conditionalP) const;

    void updateNodes(); // re-assign nodes to susceptible/infectious vectors

    double getOverestimatedP(const Cell &otherCell) const;

    std::vector<std::string> getSerotypesPresent() const;

    void addSerotypePresent(const std::string &);

    void reassignNodes(std::vector<Node> &newNodes); // re-calculate which Nodes are within cell for copy constructor

    void updateSerotypeNodeMap();

    bool operator==(const Cell &rhs) const;

    bool operator!=(const Cell &rhs) const;
};


#endif //MODEL_SINGLESEROTYPE_CELL_H
