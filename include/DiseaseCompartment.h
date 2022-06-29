//
// Created by Glen on 08/05/2020.
//

#ifndef MODEL_FMD_DISEASECOMPARTMENT_H
#define MODEL_FMD_DISEASECOMPARTMENT_H

#include <vector>
#include <random>
#include <variant>

/* The purpose of this class is as a self-contained disease compartment,
 * which can have an arbitrary number of sub-compartments and
 * supports the needed operations that need to be done with/to it.
 *
 * Each Compartment should have a Disease Status attached to it, and
 * the Status of the Compartment it flows in to. It should also record
 * the number of sub-compartments it has, and the total number of animals in it.*/


enum class Status {
    MAT, SUS, EXP, INF, REC, SUSVAC, RECVAC, CAR
};

struct CompartmentSwitch {
    Status compartmentA;
    Status compartmentB;
    double proportionA;

    CompartmentSwitch(Status a, Status b, double propA) : compartmentA(a), compartmentB(b), proportionA(propA) {}
};

struct TauLeapResult {
    std::variant<Status, CompartmentSwitch> targetCompartments;
    int numTransferred{0};

    TauLeapResult(const std::variant<Status, CompartmentSwitch> &target, const int num) : targetCompartments(target),
                                                                                          numTransferred(num) {}
};

class DiseaseCompartment {
private:
    Status m_name;
    std::variant<Status, CompartmentSwitch> m_nextCompartment;
    int m_totalPopulation;
    size_t m_numSubcompartments; // min 1!
    // assume that population goes into m_subcompartments[0]
    std::vector<int> m_subcompartments;
    std::vector<double> m_carrierAdjustedRates;

    std::vector<int> discreteSampleSubcompartments(std::mt19937 &prng, int sampleSize);

    std::vector<int> approxSampleSubcompartments(std::mt19937 &prng, int sampleSize);

public:
    DiseaseCompartment(const Status name, const std::variant<Status, CompartmentSwitch> nextCompartment,
                       const int population,
                       const int numSubcompartments);

    DiseaseCompartment(const Status name, const std::variant<Status, CompartmentSwitch> nextCompartment,
                       const int population,
                       const std::vector<int> &subcompartments);

    DiseaseCompartment(const Status name, const std::variant<Status, CompartmentSwitch> nextCompartment,
                                           const std::vector<int> &subcompartments);

    DiseaseCompartment operator+(const DiseaseCompartment &rhs);

    DiseaseCompartment operator+=(const DiseaseCompartment &rhs);

    int getTotalPopulation() const;

    int getSubPopulation(size_t subcompartmentIndex) const;

    Status getName() const;

    std::variant<Status, CompartmentSwitch> getNextCompartment() const;

    size_t getNumSubcompartments() const;

    std::vector<int> getRandomSubset(std::mt19937 &prng, int sampleSize);

    TauLeapResult tauLeap(std::mt19937 &prng, double rate);

    void addPopulation(int);

    int randomDeaths(std::mt19937 &prng, int deaths);

    void zeroCompartment();

    int orderedDeaths(int);

    void removePopulation(const std::vector<int> &subcompartmentDeaths);

    int extractFromSingleSubcompartment(size_t subcompartment, int toExtract);

    void calculateAdjustedCarrierRates(double standardTransmission);

    double getAdjustedCarrierInfectivity(double standardTransmission);

    std::vector<int> getSubCompartments() const;
};


#endif //MODEL_FMD_DISEASECOMPARTMENT_H
