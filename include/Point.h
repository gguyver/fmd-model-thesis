//
// Created by glen on 03/02/2020.
//

#ifndef MODEL_SINGLESEROTYPE_POINT_H
#define MODEL_SINGLESEROTYPE_POINT_H


class Point {
private:
    double m_xCoord{0.0};
    double m_yCoord{0.0};
public:
    explicit Point(double x = 0.0, double y = 0.0) : m_xCoord(x), m_yCoord(y) {};

    Point(int x, int y) {
        m_xCoord = static_cast<double>(x);
        m_yCoord = static_cast<double>(y);
    };

    double getXCoord() const;

    void setXCoord(double xCoord);

    double getYCoord() const;

    void setYCoord(double yCoord);

    bool operator==(const Point &rhs) const;

    bool operator!=(const Point &rhs) const;
};


#endif //MODEL_SINGLESEROTYPE_POINT_H
