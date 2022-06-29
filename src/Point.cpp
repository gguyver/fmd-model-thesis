//
// Created by glen on 03/02/2020.
//

#include <algorithm>
#include "Point.h"

bool Point::operator==(const Point &rhs) const {
    // account for float precision
    return (std::abs(m_xCoord - rhs.m_xCoord) <=
            0.00001 * std::max({1.0, std::abs(m_xCoord), std::abs(rhs.m_xCoord)}) &&
            std::abs(m_yCoord - rhs.m_yCoord) <= 0.00001 * std::max({1.0, std::abs(m_yCoord), std::abs(rhs.m_yCoord)}));
}

bool Point::operator!=(const Point &rhs) const {
    return !(rhs == *this);
}

double Point::getXCoord() const {
    return m_xCoord;
}

void Point::setXCoord(double xCoord) {
    m_xCoord = xCoord;
}

double Point::getYCoord() const {
    return m_yCoord;
}

void Point::setYCoord(double yCoord) {
    m_yCoord = yCoord;
}
