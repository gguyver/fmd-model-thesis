#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <Common.h>
#include <algorithm>
#include "Point.h"


std::vector<std::string> common::split(std::stringstream &inputString, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    while (getline(inputString, token, delim)) {
        tokens.emplace_back(token);
    }
    return (tokens);
}


/// Based on algorithm described at http://www.johndcook.com/blog/cpp_expm1/
/// \param[in]	x	Exponent values
double common::oneMinusExp(double x) {
    if (x == 0.0) {
        return 0.0;
    } else if (std::abs(x) < 1e-5) {
        return -(x + 0.5 * x * x);
    } else {
        return -(exp(x) - 1.0);
    }
}

double common::euclideanDistance(const Point &pointA, const Point &pointB) {
    double vertical{pointA.getXCoord() - pointB.getXCoord()};
    double horizontal{pointA.getYCoord() - pointB.getYCoord()};
    return std::sqrt(std::pow(vertical, 2) + std::pow(horizontal, 2));
}

bool common::checkFilePathWorks(std::string &filePath) {
    std::fstream file{filePath.c_str()};
    return file.good();
}

double common::infectionProbabilty(double transmission, double susceptibility, double kernel) {
    return (1.0 - std::exp((-transmission * susceptibility * kernel)));
}

std::string common::trimEnds(std::string const &str) {
    const auto begin = str.find_first_not_of(" \t\n\r\f\v");
    const auto end = str.find_last_not_of(" \t\n\r\f\v");
    return (str.substr(begin, begin + (end - begin)));
}

bool common::trueFalseToBool(std::string str) {
    // https://stackoverflow.com/questions/3613284/c-stdstring-to-boolean
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}
