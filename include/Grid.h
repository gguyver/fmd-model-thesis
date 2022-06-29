//
// Created by glen on 03/02/2020.
//

#ifndef MODEL_SINGLESEROTYPE_GRID_H
#define MODEL_SINGLESEROTYPE_GRID_H

// FORWARD DECLARATIONS
class Config;

struct ShipmentRecord;

struct SeedEvent;

class Reporting;

#include <unordered_map>
#include <utility>
#include <queue>
#include "Node.h"
#include "Cell.h"
#include "Matrix.h"

struct ValidatedShipment {
    ShipmentRecord *const record;
    Node *sourcePtr;
    Node *targetPtr;

    ValidatedShipment(ShipmentRecord *const rec, Node *source, Node *target) : record(
            rec), sourcePtr(source), targetPtr(target) {}
};

class Grid {
private:
    friend class Cell;

    friend class Node;

    const Config *m_config;
    Reporting *m_reporting;
    std::mt19937 *m_pRNG;
    int m_currentDay;

    static inline std::unordered_map<std::string, std::vector<std::string>> s_vaccinationRadiusIDs;
    static inline std::unordered_map<std::string, std::vector<std::string>> s_banRadiusIDs;

    double m_verticalLength{0.0}, m_horizontalLength{0.0}, m_longestSide{0.0};
    double m_longestDistance{0.0}; // distance between two opposite corners of Grid
    size_t m_numCells{0};
    size_t m_optimalNodesPerCell{0};

    std::vector<double> m_pointsX, m_pointsY, m_boundsX, m_boundsY; // minimum then maximum boundaries of X and Y planes
    Point m_lowerLeft, m_upperLeft, m_upperRight, m_lowerRight;
    // NODES
    std::vector<Node> m_allNodes;
    std::vector<Node *> m_infectiousNodes, m_detectedNodes;
    std::queue<Node *> m_vaccinationQueue;
    std::vector<Cell> m_allCells;
    std::unordered_map<std::string, Node *> m_nodeIDLookup; // gives access to Nodes via their ID
    Matrix<unsigned int> m_matrixDistanceCellToCell{0, 1, 1};
    Matrix<double> m_matrixOverestimatedProbabilities{0.0, 1,
                                                      1}; // precalculated overestimated cell-cell infection probabilities, using most susceptible serotype

    void readNodeFile();

    void readCleanNodeFile();

    void readSavedNodeFile();

    void calculateOptimalNodesPerCell();

    double calculateKernelCalls(size_t CellsOnOneSide);

    static double
    optimisationDistanceFunction(const size_t &xA, const size_t &yA, const size_t &xB, const size_t &yB,
                                 double longestSide,
                                 size_t cellsOneSide);

    void calculateDimensions();

    void createDynamicCells();

    void calculateShortestDistances();

    void calculateKernelValues();

    double calculateOverestimatedP(const Cell &cellInf, const Cell &cellSus);

    void incrementDay(int day, const std::string &runName);

    void seedInfection(const SeedEvent &seed, const std::string &runName);

    void incrementNodes();

    // get shipments that occur on specified day between unbanned nodes
    std::vector<ValidatedShipment> validateShipments(int day);

    void runShipments(const std::string &runName, int day);

    void detectedNodes(); // detects whether infected regardless of serotype

    void infectedNodes(); // updates m_infectiousNodes, should be run before sumSerotypeStatistics/Prevalence

    void vaccinationQueueNodes(); // for use with saved node state

    std::pair<int, int> runVaccination();

    void identifyReactiveVaccinationTargets();

    std::pair<int, int> runMovementBans();

    void updateMovementBans();

    void identifyMassVaccinationTargets();

    void infectionForcing(const std::string &runName, const std::string &dayString);

    void updateCells();

    void runLocalSpread(const std::string &runName, const std::string &dayName);

    double powerLawKernel(double distance);

    void initNodeRadii();

    int seedFixedSeeds(const SeedEvent &seed, std::string const &runName) const;

    int seedRandomSeeds(const SeedEvent &seed, std::string const &runName);

public:


    /// Constructor
    explicit Grid(const Config *config);

    /// Destructor
    ~Grid() = default;

    /// Copy Constructor
    Grid(const Grid &rhs); // Copy constructor

    /// Copy Assignment
    Grid &operator=(const Grid &);

    /// Move Constructor
    Grid(Grid &&other) = delete;

    /// Move Assignment
    Grid &operator=(Grid &&other) = delete;

    // run the model from start_day -> end_day
    void run(const std::string &name, bool isBurnIn);

    void saveNodeState();

};


#endif //MODEL_SINGLESEROTYPE_GRID_H
