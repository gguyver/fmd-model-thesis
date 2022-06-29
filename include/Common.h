// header for the model

#ifndef MODEL_H
#define MODEL_H

// FORWARD DECLARATIONS
class Point;

#include <map>
#include <string>
#include <vector>
#include <random>


namespace common {
    bool checkFilePathWorks(std::string &filePath);

    std::vector<std::string> split(std::stringstream &inputString, char delim);

    double oneMinusExp(double x);

    // template functions need the implementation available when linking, so in header files
    template<typename T>
    double average(const std::vector<T> &elements) {
        double sum = 0.0;
        for (T element : elements) {
            sum += double(element);
        }
        return sum / double(elements.size());
    }

    template<typename T>
    double variance_from(const std::vector<T> &elements, double val) {
        std::vector<double> sq_diff;
        sq_diff.reserve(elements.size());
        for (T element : elements) {
            double diff = double(element) - val;
            sq_diff.push_back(diff * diff);
        }
        return common::average(sq_diff);
    }

    double euclideanDistance(const Point &pointA, const Point &pointB);

    double infectionProbabilty(double transmission, double susceptibility, double kernel);

    template<typename K, typename V>
    std::vector<V> mapValuesToVector(std::map<K, V> map) {
        std::vector<V> values;
        for (const auto&[key, value] : map) {
            values.push_back(value);
        }
        return values;
    }

    template<typename T>
    T weightedSampleVector(std::mt19937 &RNG, const std::vector<T> &sampleVector, const std::vector<int> &weights) {
        std::discrete_distribution<int> dist(weights.begin(), weights.end());
        return sampleVector[dist(RNG)];
    }

    std::string trimEnds(std::string const &str);

    bool trueFalseToBool(std::string str);
}



#endif
