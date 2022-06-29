#ifndef EPI_H
#define EPI_H

// FORWARD DECLARATIONS
class Grid;

class Cell;

class Config;

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <random>
#include <map>
#include <memory>
#include "yaml.h"
#include "Point.h"
#include "SerotypeCompartmentalModel.h"

enum class Status;

struct Shipment {
    int size{0};
    bool sourceIsInfected{false};
    bool shipmentIsInfected{false};
    std::string infectionSerotype{"none"};
    std::vector<SerotypeCompartmentalModel> serotypeCompartments;
};

// Communicate result of accepting shipment back to Grid::runShipments().
struct ShipmentResult {
    bool targetInfectedPrior{false};
    bool targetNowInfected{false};
    bool fomiteTransmissionOccured{false};

};


class Node {
private:
    const std::string m_epiID{"0"};
    const Point m_location{0.0, 0.0};
    const Config *m_config;
    Grid *m_parentGrid;
    Cell *m_parentCell;
    std::mt19937 *m_pRNG;

    // Compartments
    int m_numCattle{0};
    std::vector<SerotypeCompartmentalModel> m_serotypeCompartments; // at expected small n vector is faster than map

    // Infection
    bool m_isInfected{false};
    std::string m_serotypeInfectedWith{"none"};
    std::pair<double, double> m_transmissionPars{0.0,
                                                 0.0}; // cattle-level, infectiousness, susceptibility 0 if un-infected, serotype-specific value otherwise
    bool m_infectionDetectable{true};
    bool m_hasBeenDetected{false};
    int m_daysSinceInfected{0};
    // movements
    bool m_isSentinelNode{false};
    bool m_allowedMovements{true};
    int m_daysSinceMovementBan{0};
    std::vector<std::string> m_banRadiusIDs;
    // Vaccination
    bool m_isQueuedForVaccination{false};
    bool m_hasBeenVaccinated{false};
    int m_daysSinceVaccination{0};
    std::vector<std::string> m_vaccinationRadiusIDs;

    std::vector<SerotypeCompartmentalModel>::iterator findSerotype(std::string_view serotype);

    std::vector<SerotypeCompartmentalModel>::const_iterator findSerotype(std::string_view serotype) const;

    std::pair<int, int> calculateBirthsAndDeaths() const;

    void checkIsInfected();

    void updateDSI();

    void updateDaysSinceMovementBan();

    void updateDaysSinceVaccinated();

    void updateTransmission();

public:
    Node(std::string ID, Point location, int numCattle, Grid *gridPtr, const Config *config, std::mt19937 *rng);

    Node(YAML::Node const & yaml, Grid *gridPtr, const Config *config, std::mt19937 *rng);

    [[nodiscard]] std::string getID() const;

    Point getLocation() const;

    int getNumCattle() const;

    bool isInfected() const;

    bool isInfectiousWith(const std::string &) const;

    const std::string &getSerotypeInfectedWith() const;

    bool detectInfected();

    bool getHasBeenDetected() const;

    int getDaysSinceInfected() const;

    bool getAllowedToMove() const;

    int getDaysSinceMovementBan() const;

    bool hasBeenVaccinated() const;

    int getDaysSinceVaccinated() const;

    void setIsQueuedForVaccination(bool mIsQueuedForVaccination);

    bool getIsQueuedForVaccination() const;

    void setGridIn(Grid *gridPtr);

    Cell *getCellIn() const;

    void setCellIn(Cell *cell);

    bool canBeInfectedWithSerotype(const std::string &infectingSerotype) const;

    void setMovement(bool newMovementPermission);

    void update();

    std::pair<bool, std::string>
    attemptToInfect(std::map<std::string, double> &serotypeInfectionProbabilities);

    void uninfect();

    void tauLeap();

    Shipment generateShipment(int numShipped);

    ShipmentResult acceptShipment(Shipment &incomingShipment);

    int vaccinate();

    double getCattleInfectionProb(const Node *otherNode);

    void setIsSentinel(bool);

    std::vector<int> getSumCompartments(std::string_view serotype) const;

    void readState(std::vector<std::string> const &tokens);

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Node &node);

};

#endif