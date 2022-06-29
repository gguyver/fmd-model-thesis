//
// Created by Glen on 08/05/2020.
//

#include <stdexcept>
#include <random>
#include "DiseaseCompartment.h"

// Constructor
DiseaseCompartment::DiseaseCompartment(const Status name, const std::variant<Status, CompartmentSwitch> nextCompartment,
                                       const int population, const int numSubcompartments) :
        m_name(name), m_nextCompartment(nextCompartment), m_totalPopulation(population) {
    if (std::holds_alternative<Status>(m_nextCompartment)) {
        if (m_name == std::get<Status>(m_nextCompartment)) {
            throw std::invalid_argument(
                    "DiseaseCompartment: m_nextCompartment cannot hold the same Status(es) as m_name.");
        }
    } else {
        auto *cSwitch = &std::get<CompartmentSwitch>(m_nextCompartment);
        if (m_name == cSwitch->compartmentA || m_name == cSwitch->compartmentB) {
            throw std::invalid_argument(
                    "DiseaseCompartment: m_nextCompartment cannot hold the same Status(es) as m_name.");
        }
    }

    if (numSubcompartments < 1) {
        throw std::invalid_argument("DiseaseCompartment: number of subcompartments cannot be < 1.");
    }
    if (m_totalPopulation < 0) {
        throw std::invalid_argument("DiseaseCompartment: population cannot be negative");
    }
    // assign population to first sub-compartment
    m_numSubcompartments = static_cast<size_t>(numSubcompartments);
    m_subcompartments = std::vector<int>(m_numSubcompartments, 0);
    m_subcompartments[0] = m_totalPopulation;
}

DiseaseCompartment::DiseaseCompartment(const Status name, const std::variant<Status, CompartmentSwitch> nextCompartment,
                                       const int population,
                                       const std::vector<int> &subcompartments) :
        m_name(name), m_nextCompartment(nextCompartment), m_totalPopulation(population),
        m_subcompartments(subcompartments) {
    if (std::holds_alternative<Status>(m_nextCompartment)) {
        if (m_name == std::get<Status>(m_nextCompartment)) {
            throw std::invalid_argument(
                    "DiseaseCompartment: m_nextCompartment cannot hold the same Status(es) as m_name.");
        }
    } else {
        auto *cSwitch = &std::get<CompartmentSwitch>(m_nextCompartment);
        if (m_name == cSwitch->compartmentA || m_name == cSwitch->compartmentB) {
            throw std::invalid_argument(
                    "DiseaseCompartment: m_nextCompartment cannot hold the same Status(es) as m_name.");
        }
    }
    if (m_totalPopulation < 0) {
        throw std::invalid_argument("DiseaseCompartment: population cannot be negative");
    }
    m_numSubcompartments = subcompartments.size();
}

DiseaseCompartment::DiseaseCompartment(const Status name, const std::variant<Status, CompartmentSwitch> nextCompartment,
                                       const std::vector<int> &subcompartments) :
        m_name(name), m_nextCompartment(nextCompartment),
        m_subcompartments(subcompartments) {
    if (std::holds_alternative<Status>(m_nextCompartment)) {
        if (m_name == std::get<Status>(m_nextCompartment)) {
            throw std::invalid_argument(
                    "DiseaseCompartment: m_nextCompartment cannot hold the same Status(es) as m_name.");
        }
    } else {
        auto *cSwitch = &std::get<CompartmentSwitch>(m_nextCompartment);
        if (m_name == cSwitch->compartmentA || m_name == cSwitch->compartmentB) {
            throw std::invalid_argument(
                    "DiseaseCompartment: m_nextCompartment cannot hold the same Status(es) as m_name.");
        }
    }
    m_totalPopulation = 0;
    for(const auto& sub : m_subcompartments) {
        m_totalPopulation += sub;
    }
    if (m_totalPopulation < 0) {
        throw std::invalid_argument("DiseaseCompartment: population cannot be negative");
    }
    m_numSubcompartments = subcompartments.size();
}

int DiseaseCompartment::getTotalPopulation() const {
    return m_totalPopulation;
}

int DiseaseCompartment::getSubPopulation(size_t subcompartmentIndex) const {
    // range check
    if (subcompartmentIndex > (m_numSubcompartments - 1)) {
        throw std::out_of_range("DiseaseCompartment subcompartment does not exist");
    }
    return m_subcompartments[subcompartmentIndex];
}

Status DiseaseCompartment::getName() const {
    return m_name;
}

std::variant<Status, CompartmentSwitch> DiseaseCompartment::getNextCompartment() const {
    return m_nextCompartment;
}

size_t DiseaseCompartment::getNumSubcompartments() const {
    return m_numSubcompartments;
}

void DiseaseCompartment::addPopulation(int newPopulation) {
    if (newPopulation < 0) {
        throw std::invalid_argument("DiseaseCompartment: Do not use addPopulation with negative integers");
    } else if (newPopulation > 0) {
        m_totalPopulation += newPopulation;
        m_subcompartments[0] += newPopulation;
    }
}

int DiseaseCompartment::randomDeaths(std::mt19937 &prng, int deaths) {
    if (deaths < 0) {
        throw std::invalid_argument("DiseaseCompartment::randomDeaths cannot take a negative argument");
    }
    if (deaths >= m_totalPopulation) {
        deaths -= m_totalPopulation;
        zeroCompartment();
        return deaths;
    } else {
        if (m_numSubcompartments > 1) {
            getRandomSubset(prng, deaths);
        } else {
            m_subcompartments[0] -= deaths;
            m_totalPopulation -= deaths;
        }
        return 0;
    }
}

int DiseaseCompartment::orderedDeaths(int deaths) {
    // intended for use with
    if (deaths < 0) {
        throw std::invalid_argument("DiseaseCompartment::orderedDeaths cannot take a negative argument");
    }
    if (deaths >= m_totalPopulation) {
        // return excess deaths
        deaths -= m_totalPopulation;
        zeroCompartment();
        return deaths;
    } else {
        m_totalPopulation -= deaths;
        if (m_numSubcompartments > 1) {
            // reverse order
            for (auto rit = m_subcompartments.rbegin(); rit != m_subcompartments.rend() && deaths > 0; rit++) {
                int remainder = std::max(0, (*rit) - deaths);
                deaths = std::max(0, deaths - (*rit));
                (*rit) = remainder;
            }
        } else {
            m_subcompartments[0] -= deaths;
        }

        return 0;
    }
}

TauLeapResult DiseaseCompartment::tauLeap(std::mt19937 &prng, double rate) {
    // Tau-leap algorithm to approximate gillespie algorithm
    if (rate < 0.0) {
        throw std::invalid_argument("DiseaseCompartment::tauLeap rate must be > 0.0");
    }

    /* rate multiplied by number of subcompartments (n) so travel through each compartment
         at N-times the speed to spend the same average time in the overall compartment */
    std::vector<double> subcompartmentRates(m_numSubcompartments, 0.0);
    for (size_t i = 0; i < m_numSubcompartments; ++i) {
        subcompartmentRates[i] = m_subcompartments[i] * rate * static_cast<double>(m_numSubcompartments);
    }

    /* Draw number of events from poisson distribution, underlying this is an exponential distribution.
     * As m_numSubcompartments increases, the overall distribution becomes a gamma distribution */
    int toNextCompartment{0};
    for (size_t i = 0; i < m_numSubcompartments; ++i) {
        double r{subcompartmentRates[i]};
        if (r > 0.0) {
            std::poisson_distribution<int> pois(r);
            // make sure events doesn't exceed subcompartment population
            int events{std::min(m_subcompartments[i], pois(prng))};

            /* move population to next subcompartment, except if last subcompartment
             * where actually return that value with reference to next DiseaseCompartment
             * to add it to */
            m_subcompartments[i] -= events;
            if (i < (m_numSubcompartments - 1)) {
                m_subcompartments[(i + 1)] += events;
            } else {
                toNextCompartment += events;
                m_totalPopulation -= events;
            }
        }
    }
    return TauLeapResult(m_nextCompartment, toNextCompartment);
}

void DiseaseCompartment::zeroCompartment() {
    m_totalPopulation = 0;
    std::fill(m_subcompartments.begin(), m_subcompartments.end(), 0);
}

std::vector<int> DiseaseCompartment::getRandomSubset(std::mt19937 &prng, int sampleSize) {
    /* Need to sample without replacement. Unfortunately there is no easy
     * way to do this in standard C++. The two options I can see are:
     * - Create vector of length population, with values equating to indices of the subcompartments
     *      and use std::sample to sample this without replacement
     * - Use a for loop with std::discrete_distribution, although
     *      this involves creating the distribution every loop */
    if (sampleSize < 0) {
        throw std::invalid_argument("DiseaseCompartment::getRandomSubset cannot sample a negative size sample");
    }

    std::vector<int> subset(m_numSubcompartments, 0);
    if (sampleSize < m_totalPopulation) {
        if (sampleSize <= 5) {
            subset = discreteSampleSubcompartments(prng, sampleSize);
        } else {
            subset = approxSampleSubcompartments(prng, sampleSize);
        }
    } else {
        // setting to m_totalPopulation will not be different to m_subcompartments so just return that
        subset = m_subcompartments;
        zeroCompartment();
    }
    return subset;
}

DiseaseCompartment DiseaseCompartment::operator+(const DiseaseCompartment &rhs) {
    m_totalPopulation += rhs.m_totalPopulation;
    for (size_t i = 0; i < m_subcompartments.size(); ++i) {
        m_subcompartments[i] += rhs.m_subcompartments[i];
    }
    return *this;
}

DiseaseCompartment DiseaseCompartment::operator+=(const DiseaseCompartment &rhs) {
    return (*this + rhs);
}

int DiseaseCompartment::extractFromSingleSubcompartment(size_t subcompartment, int toExtract) {
    /* This function is basically for the rare times when there is no need to randomise the extract,
     * as you know exactly where the population should come from. E.g. when taking from Status::SUS
     * which has no need for multiple subcompartments */
    if (subcompartment > m_numSubcompartments - 1) {
        throw std::invalid_argument("Cannot extract from subcompartment that does not exist");
    } else if (toExtract < 0 || toExtract > m_subcompartments[subcompartment]) {
        throw std::invalid_argument("toExtract is an invalid value.");
    }
    if (toExtract > 0) {
        m_subcompartments[subcompartment] -= toExtract;
        m_totalPopulation -= toExtract;
        return toExtract;
    } else {
        return 0;
    }
}

void DiseaseCompartment::removePopulation(const std::vector<int> &subcompartmentDeaths) {
    if (subcompartmentDeaths.size() != m_numSubcompartments) {
        throw std::invalid_argument(
                "DiseaseCompartment::removePopulation() : provided vector is not equal in length to subcompartments.");
    }
    for (size_t i = 0; i < m_subcompartments.size(); ++i) {
        if (subcompartmentDeaths[i] > m_subcompartments[i]) {
            throw std::invalid_argument(
                    "DiseaseCompartment::deaths() : subcompartment removePopulation greater than subcompartment population");
        }
        m_subcompartments[i] -= subcompartmentDeaths[i];
        m_totalPopulation -= subcompartmentDeaths[i];
    }
}

std::vector<int> DiseaseCompartment::discreteSampleSubcompartments(std::mt19937 &prng, int sampleSize) {
    if (sampleSize > m_totalPopulation) {
        throw std::invalid_argument("discrete sample sampleSize cannot be > population");
    }
    std::vector<int> sampleCounts(m_numSubcompartments, 0);
    if (m_numSubcompartments > 1) {
        for (int i = 0; i < sampleSize; ++i) {
            std::discrete_distribution<size_t> dist(m_subcompartments.begin(), m_subcompartments.end());
            size_t selected{dist(prng)};
            ++sampleCounts[selected];
            --m_subcompartments[selected];
            --m_totalPopulation;
        }
    } else {
        sampleCounts[0] += sampleSize;
        m_subcompartments[0] -= sampleSize;
        m_totalPopulation -= sampleSize;
    }
    return sampleCounts;
}

std::vector<int> DiseaseCompartment::approxSampleSubcompartments(std::mt19937 &prng, int sampleSize) {
    if (sampleSize > m_totalPopulation) {
        throw std::invalid_argument("approx sample sampleSize cannot be > population");
    }
    std::vector<int> sampleCounts(m_numSubcompartments, 0);
    if (m_numSubcompartments > 1) {
        double sampleProportion{double(sampleSize) / double(m_totalPopulation)};
        int sumSampled{0};
        for (size_t i = 0; i < m_subcompartments.size() && sumSampled < sampleSize; ++i) {
            int thisSample{int(std::floor(sampleProportion * double(m_subcompartments[i])))};
            if (thisSample > 0) {
                m_subcompartments[i] -= thisSample;
                m_totalPopulation -= thisSample;
                sampleCounts[i] += thisSample;
                sumSampled += thisSample;
            }

        }
        if ((sampleSize - sumSampled) > 0) {
            std::vector<int> remainder{discreteSampleSubcompartments(prng, (sampleSize - sumSampled))};
            for (size_t i = 0; i < sampleCounts.size(); ++i) {
                sampleCounts[i] += remainder[i];
            }
        }
    } else {
        sampleCounts[0] += sampleSize;
        m_subcompartments[0] -= sampleSize;
        m_totalPopulation -= sampleSize;
    }
    return sampleCounts;
}

void DiseaseCompartment::calculateAdjustedCarrierRates(double standardTransmission) {
    if (m_carrierAdjustedRates.empty()) {
        // count number of 'units' of infectivity per subcompartment
        size_t subCompartmentDivider{0};
        for (size_t i = 1; i <= m_numSubcompartments; ++i) {
            subCompartmentDivider += i;
        }
        double smallestUnit{standardTransmission / static_cast<double>(subCompartmentDivider)};
        for (size_t j = 0; j < m_subcompartments.size(); ++j) {
            double subAdjustedInfectivity{(static_cast<double>((m_numSubcompartments - j)) * smallestUnit)};
            m_carrierAdjustedRates.push_back(subAdjustedInfectivity);
        }
    }
}

double DiseaseCompartment::getAdjustedCarrierInfectivity(double standardTransmission) {
    calculateAdjustedCarrierRates(standardTransmission);
    double adjustedInfectivity{0.0};
    for (size_t i = 0; i < m_subcompartments.size(); ++i) {
        adjustedInfectivity += m_subcompartments[i] * m_carrierAdjustedRates[i];
    }
    return adjustedInfectivity;
}

std::vector<int> DiseaseCompartment::getSubCompartments() const {
    return m_subcompartments;
}


