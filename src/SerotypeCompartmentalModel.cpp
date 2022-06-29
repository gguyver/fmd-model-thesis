//
// Created by Glen on 12/05/2020.
//

#include <stdexcept>
#include <algorithm>
#include "SerotypeCompartmentalModel.h"
#include "Config.h"
#include "Grid.h"

SerotypeCompartmentalModel::SerotypeCompartmentalModel(const Config *config, const std::string &serotype, int pop)
        : m_config{config}, m_serotype{serotype} {
    if (pop < 0) {
        throw std::invalid_argument("SerotypeCompartmentalModel cannot be initialised with a negative population");
    }
    // Hard coding this because it won't be changing at run-time.

    m_population = pop;
    m_compartments.emplace_back(Status::MAT, Status::SUS, 0, 3);
    m_compartments.emplace_back(Status::SUS, Status::EXP, pop, 1);
    m_compartments.emplace_back(Status::EXP, Status::INF, 0, 1);
    m_compartments.emplace_back(Status::INF, CompartmentSwitch(Status::CAR, Status::REC,
                                                               m_config->m_serotypeInfectionParameters.at(
                                                                       serotype).proportionCarriers), 0, 1);
    m_compartments.emplace_back(Status::REC, Status::SUS, 0, 20);
    m_compartments.emplace_back(Status::SUSVAC, Status::SUS, 0, 20);
    m_compartments.emplace_back(Status::RECVAC, Status::SUS, 0, 20);
    m_compartments.emplace_back(Status::CAR, Status::REC, 0, 5);

    for (const auto &comp : m_compartments) {
        m_numSubcompartments += comp.getNumSubcompartments();
    }

    m_infectedStatuses = {Status::EXP, Status::INF};
    m_infectiousStatuses = {Status::INF};
    m_immuneStatuses = {Status::REC, Status::SUSVAC, Status::RECVAC};
    m_asymptomaticAdultStatuses = {Status::SUS, Status::EXP, Status::REC, Status::SUSVAC, Status::RECVAC};
}

SerotypeCompartmentalModel::SerotypeCompartmentalModel(const Config *config, const std::string &serotype, const YAML::Node& yaml)
        : m_config{config}, m_serotype{serotype} {

    m_compartments.emplace_back(Status::MAT, Status::SUS, yaml["MAT"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::SUS, Status::EXP, yaml["SUS"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::EXP, Status::INF, yaml["EXP"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::INF, CompartmentSwitch(Status::CAR, Status::REC,
                                                               m_config->m_serotypeInfectionParameters.at(
                                                                       serotype).proportionCarriers), yaml["INF"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::REC, Status::SUS, yaml["REC"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::SUSVAC, Status::SUS, yaml["SUSVAC"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::RECVAC, Status::SUS, yaml["RECVAC"].as<std::vector<int>>());
    m_compartments.emplace_back(Status::CAR, Status::REC, yaml["CAR"].as<std::vector<int>>());

    m_population = 0;
    for (const auto &comp : m_compartments) {
        m_population += comp.getTotalPopulation();
        m_numSubcompartments += comp.getNumSubcompartments();
    }
    if (m_population < 0) {
        throw std::invalid_argument("SerotypeCompartmentalModel cannot be initialised with a negative population");
    }

    m_infectedStatuses = {Status::EXP, Status::INF};
    m_infectiousStatuses = {Status::INF};
    m_immuneStatuses = {Status::REC, Status::SUSVAC, Status::RECVAC};
    m_asymptomaticAdultStatuses = {Status::SUS, Status::EXP, Status::REC, Status::SUSVAC, Status::RECVAC};
}

SerotypeCompartmentalModel::SerotypeCompartmentalModel(const Config *pConfig, const std::string &serotype,
                                                       int population,
                                                       const std::vector<DiseaseCompartment> &compartments) :
        m_config(pConfig), m_serotype(serotype), m_population(population),
        m_compartments(compartments) {

    for (const auto &comp : m_compartments) {
        m_numSubcompartments += comp.getNumSubcompartments();
    }

    m_infectedStatuses = {Status::EXP, Status::INF};
    m_infectiousStatuses = {Status::INF};
    m_immuneStatuses = {Status::REC, Status::SUSVAC, Status::RECVAC};
    m_asymptomaticAdultStatuses = {Status::SUS, Status::EXP, Status::REC, Status::SUSVAC, Status::RECVAC};
}

std::vector<DiseaseCompartment>::iterator SerotypeCompartmentalModel::findCompartment(const Status &status) {
    auto it = std::find_if(m_compartments.begin(),
                           m_compartments.end(),
                           [&](DiseaseCompartment &comp) { return comp.getName() == status; });
    if (it != m_compartments.end()) {
        return it;
    } else {
        // should not be able to search for a nonexistent Status
        throw std::invalid_argument("SerotypeCompartmentModel cannot find non-existent Status");
    }
}

std::vector<DiseaseCompartment>::const_iterator
SerotypeCompartmentalModel::findCompartment(const Status &status) const {
    auto it = std::find_if(m_compartments.begin(),
                           m_compartments.end(),
                           [&](const DiseaseCompartment &comp) { return comp.getName() == status; });
    if (it != m_compartments.end()) {
        return it;
    } else {
        // should not be able to search for a nonexistent Status
        throw std::invalid_argument("SerotypeCompartmentModel cannot find non-existent Status");
    }
}

int SerotypeCompartmentalModel::getTotalPopulation() const {
    return m_population;
}

int SerotypeCompartmentalModel::getCompartmentPopulation(const Status &status) const {
    auto compartment = findCompartment(status);
    return (*compartment).getTotalPopulation();
}

bool SerotypeCompartmentalModel::getIsInfected() const {
    for (const auto &status : m_infectedStatuses) {
        auto it = findCompartment(status);
        if ((*it).getTotalPopulation() > 0) {
            return true;
        }
    }
    return false;
}

int SerotypeCompartmentalModel::getSumInfectious() const {
    int sum{0};
    for (const auto &status : m_infectiousStatuses) {
        auto it = findCompartment(status);
        sum += (*it).getTotalPopulation();
    }
    return sum;
}

void SerotypeCompartmentalModel::tauLeap(std::mt19937 &prng) {
    /* This function only handles non-population-changing events. Births and deaths will
     * be handled in the wrapper class, as they need to be synchronised between SCMs objects
     * so that each serotype is referring to the same total population. */

    // calculate per-capita rates
    const auto rates = calculateTauLeapRates();

    // calculate results before moving population between compartments
    const auto results = calculateTauLeapResults(prng, rates);

    // resolve population moving between compartments
    if (not results.empty()) {
        resolveTauLeapResults(prng, results);
    }

}

std::vector<TauLeapRate> SerotypeCompartmentalModel::calculateTauLeapRates() {
    // TODO: Double-check maths for CAR infection.
    auto pars{m_config->m_serotypeInfectionParameters.at(m_serotype)};
    std::vector<TauLeapRate> rates{
            {Status::MAT,    pars.maternalImmunityDecayRate},
            {Status::SUS,    (pars.beta * double(getSumInfectious()) / double(m_population)) +
                             (getAdjustedCarrierTransmission(pars.carrierTransmission) /
                              static_cast<double>(m_population))},
            {Status::EXP,    pars.symptomaticRate},
            {Status::INF,    pars.recoveryRate},
            {Status::REC,    pars.waningNatural},
            {Status::SUSVAC, pars.waningVaccine},
            {Status::RECVAC, std::min(pars.waningNatural, pars.waningVaccine)},
            {Status::CAR,    pars.waningCarriers}
    };
    return rates;
}

std::vector<TauLeapResult>
SerotypeCompartmentalModel::calculateTauLeapResults(std::mt19937 &prng,
                                                    std::vector<TauLeapRate> const &rates) {
    std::vector<TauLeapResult> results;
    for (const auto &rate : rates) {
        auto it = findCompartment(rate.compartment);
        if ((*it).getTotalPopulation() > 0) {
            // DiseaseCompartment::tauLeap() extracts population from compartment
            const auto result = (*it).tauLeap(prng, rate.rate);
            if (result.numTransferred > 0) {
                results.push_back(result);
            }
        }
    }
    return results;
}

void SerotypeCompartmentalModel::resolveTauLeapResults(std::mt19937 &prng, std::vector<TauLeapResult> const &results) {
    const auto resolveResult = [&](TauLeapResult const &res) {
        if (res.numTransferred > 0) {
            if (std::holds_alternative<Status>(res.targetCompartments)) {
                auto it = findCompartment(std::get<Status>(res.targetCompartments));
                (*it).addPopulation(res.numTransferred);
            } else {
                CompartmentSwitch const *cSwitch{&(std::get<CompartmentSwitch>(res.targetCompartments))};
                auto itA = findCompartment(cSwitch->compartmentA);
                auto itB = findCompartment(cSwitch->compartmentB);
                std::binomial_distribution<int> binomialDistribution(res.numTransferred, cSwitch->proportionA);
                int numToA{binomialDistribution(prng)};
                itA->addPopulation(numToA);
                itB->addPopulation((res.numTransferred - numToA));
            }
        }
    };
    std::for_each(results.begin(), results.end(), resolveResult);
}

void SerotypeCompartmentalModel::uninfect() {
    for (const auto &status : m_infectedStatuses) {
        auto it = findCompartment(status);
        int pop{(*it).getTotalPopulation()};
        (*it).zeroCompartment();
        auto dest_it = findCompartment(Status::REC);
        (*dest_it).addPopulation(pop);
    }
}

int SerotypeCompartmentalModel::infectionMortality(std::mt19937 &prng) {
    int sumDeaths{0};
    double deathRate{m_config->m_serotypeInfectionParameters.at(m_serotype).associatedMortalityRate};
    for (const auto &status : m_infectiousStatuses) {
        auto it = findCompartment(status);
        double pop{double((*it).getTotalPopulation())};
        if (pop > 0.0) {
            std::poisson_distribution<int> dist(deathRate * (*it).getTotalPopulation());
            int deaths{dist(prng)};
            int remainder = (*it).orderedDeaths(deaths);
            sumDeaths += (deaths - remainder);
        }
    }
    m_population -= sumDeaths;
    return sumDeaths;
}

void SerotypeCompartmentalModel::imposeDeaths(std::mt19937 &prng, int deaths) {
    /* To avoid expensive creation of random distributions, will use std::sample with a vector of (indices * population) */
    // sampleSubcompartmentsWithoutReplacement() handles deaths < 0
    std::vector<int> sampleCounts{sampleSubcompartmentsWithoutReplacement(prng, deaths)};
    removeSampleFromCompartments(sampleCounts);
}

int SerotypeCompartmentalModel::vaccinate(std::mt19937 &prng, double efficacy) {
    /* Assuming that actual efficacy has been calculated further up the chain,
     * and that this value is fixed. Also that [0.0, 1.0] */

    // Susceptible animals vaccinated
    auto susVac_it = findCompartment(Status::SUSVAC);
    auto susIt = findCompartment(Status::SUS);
    int susVaccinated = static_cast<int>(std::floor(double((*susIt).getTotalPopulation()) * efficacy));
    if (susVaccinated > 0) {
        movePopBetween(prng, susIt, susVac_it, susVaccinated);
    }

    // Recovered animals vaccinated
    auto recVac_it = findCompartment(Status::RECVAC);
    auto recIt = findCompartment(Status::REC);
    int recVaccinated = static_cast<int>(std::floor(double((*recIt).getTotalPopulation()) * efficacy));
    if (recVaccinated > 0) {
        movePopBetween(prng, recIt, recVac_it, recVaccinated);
    }

    return getSumAsymptomaticAdults();
}

void SerotypeCompartmentalModel::movePopBetween(std::mt19937 &prng, std::vector<DiseaseCompartment>::iterator &from,
                                                std::vector<DiseaseCompartment>::iterator &to, int population) {
    // TODO: don't assume population <= (*from).getTotalPopulation()
    (*from).getRandomSubset(prng, population); // deliberately discard results
    (*to).addPopulation(population);
}

SerotypeCompartmentalModel SerotypeCompartmentalModel::sampleForShipment(std::mt19937 &prng, int sampleSize) {
    // sampleSubcompartmentsWithoutReplacement() handles samples < 0 etc.
    std::vector<int> sampleCounts{sampleSubcompartmentsWithoutReplacement(prng, sampleSize)};
    removeSampleFromCompartments(sampleCounts);
    std::vector<DiseaseCompartment> sampleCompartments{sampleToCompartments(sampleCounts)};

    int actualSampleSize{std::accumulate(sampleCounts.begin(), sampleCounts.end(), 0)};
    SerotypeCompartmentalModel scm(m_config, m_serotype, actualSampleSize, sampleCompartments);
    return scm;
}

std::vector<int>
SerotypeCompartmentalModel::sampleSubcompartmentsWithoutReplacement(std::mt19937 &prng, int sampleSize) {
    if (sampleSize < 1) {
        throw std::invalid_argument("SCM::sampleSubcompartmentsWithoutReplacement() sampleSize must be >= 1");
    }
    std::vector<int> populations{getAllCompartments()};
    std::vector<int> sampleCounts;

    if (sampleSize < m_population) {
        /* discreteSample is more accurate, but unfortunately takes a long time as sampleSize increases
         * because it creates a discrete_distribution sampleSize times. approxSample is slightly less accurate
         * but much faster at larger sample sizes. */
        if (sampleSize <= 5) {
            sampleCounts = discreteSample(prng, populations, sampleSize);
        } else {
            sampleCounts = approxSample(prng, populations, sampleSize);
        }
        return sampleCounts;
    } else {
        return populations;
    }
}


void SerotypeCompartmentalModel::births(std::mt19937 &prng, int births) {
    if (births < 0) {
        throw std::invalid_argument("SCM::births() cannot accept a negative number of births");
    } else if (births > 0) {
        double sumImmune{double(getSumImmune())};
        if (m_config->m_maternalImmunityIsActive && sumImmune > 0.0) {
            double probMaternallyImmune = sumImmune / double(m_population - getSumMaternal());
            std::binomial_distribution<int> dist(births, probMaternallyImmune);
            int numMat{dist(prng)};
            addToCompartment(Status::MAT, numMat);
            addToCompartment(Status::SUS, (births - numMat));
        } else {
            addToCompartment(Status::SUS, births);
        }
        m_population += births;
    }


}

std::string SerotypeCompartmentalModel::getSerotype() const {
    return m_serotype;
}

void SerotypeCompartmentalModel::addToCompartment(const Status &status, int num) {
    auto it{findCompartment(status)};
    (*it).addPopulation(num);
}

int SerotypeCompartmentalModel::getSumImmune() {
    int sum{0};
    for (const auto &status : m_immuneStatuses) {
        auto it{findCompartment(status)};
        sum += (*it).getTotalPopulation();
    }
    return sum;
}

SerotypeCompartmentalModel
SerotypeCompartmentalModel::operator+(const SerotypeCompartmentalModel &rhs) {
    // Can rely on vector being in same order
    for (size_t i = 0; i < m_compartments.size(); ++i) {
        m_compartments[i] += rhs.m_compartments[i];
    }
    m_population += rhs.m_population;
    return *this;
}

SerotypeCompartmentalModel SerotypeCompartmentalModel::operator+=(const SerotypeCompartmentalModel &rhs) {
    return (*this + rhs);
}

int SerotypeCompartmentalModel::getSumSusceptible() const {
    auto it{findCompartment(Status::SUS)};
    return (*it).getTotalPopulation();
}

bool SerotypeCompartmentalModel::infect(int numToInfect) {
    auto sus{findCompartment(Status::SUS)};
    auto next{findCompartment(Status::EXP)};
    int infected = (*sus).extractFromSingleSubcompartment(0, numToInfect);
    if (infected > 0) {
        (*next).addPopulation(numToInfect);
        return true;
    } else {
        return false;
    }
}

std::vector<int> SerotypeCompartmentalModel::getAllCompartments() {
    // trying to avoid double allocation of vector<int>
    std::vector<int> allSubcompartments(m_numSubcompartments, 0);
    size_t traversed{0};
    for (auto &m_compartment : m_compartments) {
        for (size_t j = 0; j < m_compartment.getNumSubcompartments(); ++j) {
            allSubcompartments.at(traversed) = m_compartment.getSubPopulation(j);
            ++traversed;
        }
    }
    return allSubcompartments;
}

std::vector<int>
SerotypeCompartmentalModel::discreteSample(std::mt19937 &prng, std::vector<int> &populations, int size) {
    // WARNING: Will change vector which is passed by reference. Ok for use at the moment, as vector is copy
    std::vector<int> sampleCounts(populations.size(), 0);
    for (int i = 0; i < size; ++i) {
        std::discrete_distribution<size_t> dist(populations.begin(), populations.end());
        size_t selected{dist(prng)};
        --populations[selected];
        ++sampleCounts[selected];
    }
    return sampleCounts;
}

std::vector<int> SerotypeCompartmentalModel::approxSample(std::mt19937 &prng, std::vector<int> &populations, int size) {
    // WARNING: Will change vector which is passed by reference. Ok for use at the moment, as vector is copy
    // get proportion of each compartment, floored so know isn't going over. Then discrete method for remaining few.
    double sampleProportion{double(size) / double(m_population)};
    int sumSampled{0};
    std::vector<int> sampleCounts(populations.size(), 0);
    for (size_t i = 0; i < populations.size() && sumSampled < size; ++i) {
        int thisSample{int(std::floor(sampleProportion * double(populations[i])))};
        populations[i] -= thisSample;
        sumSampled += thisSample;
        sampleCounts[i] += thisSample;
    }
    if (sumSampled > size) {
        throw std::runtime_error("approxSample has sampled more than required. Something has gone wrong.");
    }
    if ((size - sumSampled) > 0) {
        std::vector<int> remainder{discreteSample(prng, populations, (size - sumSampled))};
        for (size_t i = 0; i < populations.size(); ++i) {
            sampleCounts[i] += remainder[i];
        }
    }
    return sampleCounts;
}

std::vector<DiseaseCompartment>
SerotypeCompartmentalModel::sampleToCompartments(const std::vector<int> &subcompartments) {
    size_t traversed{0};
    std::vector<DiseaseCompartment> compartments;
    for (const auto &comp : m_compartments) {
        std::vector<int> subvector{subcompartments.begin() + static_cast<long long int>(traversed),
                                   subcompartments.begin() + static_cast<long long int>(traversed) +
                                   static_cast<long long int>(comp.getNumSubcompartments())};
        int pop{std::accumulate(subvector.begin(), subvector.end(), 0)};
        compartments.emplace_back(comp.getName(), comp.getNextCompartment(), pop, subvector);
        traversed += comp.getNumSubcompartments();
    }
    return compartments;
}

void SerotypeCompartmentalModel::removeSampleFromCompartments(const std::vector<int> &sample) {
    size_t traversed{0};
    int sumToRemove{std::accumulate(sample.begin(), sample.end(), 0)};
    int sumRemoved{0};
    for (size_t i = 0; i < m_compartments.size() && sumRemoved < sumToRemove; ++i) {
        std::vector<int> subvector{sample.begin() + static_cast<long long int>(traversed),
                                   sample.begin() + static_cast<long long int>(traversed) +
                                   static_cast<long long int>(m_compartments[i].getNumSubcompartments())};
        int pop{std::accumulate(subvector.begin(), subvector.end(), 0)};
        if (pop > 0) {
            m_compartments[i].removePopulation(subvector);
            m_population -= pop;
            sumRemoved += pop;
        }
        traversed += m_compartments[i].getNumSubcompartments();
    }
}

int SerotypeCompartmentalModel::getSumMaternal() const {
    auto it{findCompartment(Status::MAT)};
    return (*it).getTotalPopulation();
}

std::vector<int> SerotypeCompartmentalModel::getCompartmentSums() const {
    std::vector<int> compartmentSums(m_compartments.size(), 0);
    std::transform(m_compartments.begin(), m_compartments.end(), compartmentSums.begin(),
                   [](const DiseaseCompartment &compartment) {
                       return compartment.getTotalPopulation();
                   });
    return compartmentSums;
}

int SerotypeCompartmentalModel::getSumCarriers() const {
    auto it{findCompartment(Status::CAR)};
    return (*it).getTotalPopulation();
}

double SerotypeCompartmentalModel::getAdjustedCarrierTransmission(double standardTransmission) {
    auto it{findCompartment(Status::CAR)};
    if (m_config->m_decliningCarrier) {
        return (*it).getAdjustedCarrierInfectivity(standardTransmission);
    } else {
        return (*it).getTotalPopulation() * standardTransmission;
    }
}

int SerotypeCompartmentalModel::getSumAsymptomaticAdults() {
    int sum{0};
    for (const auto &status : m_asymptomaticAdultStatuses) {
        auto it{findCompartment(status)};
        sum += (*it).getTotalPopulation();
    }
    return sum;
}

std::vector<int> SerotypeCompartmentalModel::getSubCompartments(const Status &status) const {
    auto compartment = findCompartment(status);
    return (*compartment).getSubCompartments();
}
