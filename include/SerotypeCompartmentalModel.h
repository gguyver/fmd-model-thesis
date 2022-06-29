//
// Created by Glen on 12/05/2020.
//

#ifndef MODEL_FMD_SEROTYPECOMPARTMENTALMODEL_H
#define MODEL_FMD_SEROTYPECOMPARTMENTALMODEL_H


#include "DiseaseCompartment.h"
#include "yaml.h"

// Forward declarations
class Grid;

class Config;

/* The purpose of this class is to handle the actual compartmental model stuff, which
 * will be wrapped in another class CompartmentalModel, one for each serotype */

struct TauLeapRate {
    Status compartment;
    double rate;

    TauLeapRate(Status const &status, double const rate) : compartment(status), rate(rate) {}
};

class SerotypeCompartmentalModel {
private:
    SerotypeCompartmentalModel(const Config *pConfig, const std::string &serotype, int population,
                               const std::vector<DiseaseCompartment> &compartments);

    const Config *m_config;

    const std::string m_serotype;
    int m_population{0};
    std::vector<DiseaseCompartment> m_compartments;
    size_t m_numSubcompartments{0};


    std::vector<Status> m_infectedStatuses;
    std::vector<Status> m_infectiousStatuses;
    std::vector<Status> m_immuneStatuses;
    std::vector<Status> m_asymptomaticAdultStatuses;

    [[nodiscard]] std::vector<DiseaseCompartment>::iterator findCompartment(const Status &);

    [[nodiscard]] std::vector<DiseaseCompartment>::const_iterator findCompartment(const Status &) const;

    void movePopBetween(std::mt19937 &prng, std::vector<DiseaseCompartment>::iterator &from,
                        std::vector<DiseaseCompartment>::iterator &to, int population);

    std::vector<int> sampleSubcompartmentsWithoutReplacement(std::mt19937 &prng, int sampleSize);

    void addToCompartment(const Status &, int);

    std::vector<int> discreteSample(std::mt19937 &prng, std::vector<int> &populations, int size);

    std::vector<int> approxSample(std::mt19937 &prng, std::vector<int> &populations, int size);

    std::vector<int> getAllCompartments();

    std::vector<DiseaseCompartment> sampleToCompartments(const std::vector<int> &subcompartments);

    std::vector<TauLeapRate> calculateTauLeapRates();

    std::vector<TauLeapResult>
    calculateTauLeapResults(std::mt19937 &prng, std::vector<TauLeapRate> const &rates);

    void resolveTauLeapResults(std::mt19937 &prng, std::vector<TauLeapResult> const &results);

public:
    SerotypeCompartmentalModel(const Config *config, const std::string &serotype, int pop);

    SerotypeCompartmentalModel(const Config *config, const std::string &serotype, const YAML::Node& yaml);

    SerotypeCompartmentalModel operator+(const SerotypeCompartmentalModel &rhs);

    SerotypeCompartmentalModel operator+=(const SerotypeCompartmentalModel &rhs);

    std::string getSerotype() const;

    int getTotalPopulation() const;

    int getCompartmentPopulation(const Status &) const;

    bool getIsInfected() const;

    int getSumSusceptible() const;

    int getSumInfectious() const;

    void tauLeap(std::mt19937 &prng); // only handles non-death events

    int vaccinate(std::mt19937 &prng, double efficacy);

    bool infect(int numToInfect);

    void uninfect();

    SerotypeCompartmentalModel sampleForShipment(std::mt19937 &prng, int sampleSize);

    /* Births and natural mortality should be dealt with outside the SerotypeCompartmentalModel.
     * As they are not serotype-specific and need to be synchronised to keep populations the same.
     *
     * Also serotypeMortality needs to return number dead to allow synchronisation for other SCMs*/

    int infectionMortality(std::mt19937 &prng);

    void imposeDeaths(std::mt19937 &prng, int deaths); // to synchronise SCM populations

    void births(std::mt19937 &prng, int births);

    int getSumImmune();

    void removeSampleFromCompartments(const std::vector<int> &sample);

    int getSumMaternal() const;

    std::vector<int> getCompartmentSums() const;

    int getSumCarriers() const;

    double getAdjustedCarrierTransmission(double standardTransmission);

    int getSumAsymptomaticAdults();

    std::vector<int> getSubCompartments(const Status &status) const;
};


#endif //MODEL_FMD_SEROTYPECOMPARTMENTALMODEL_H
