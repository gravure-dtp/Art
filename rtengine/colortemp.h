/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _COLORTEMP_
#define _COLORTEMP_

#include <cmath>
#include <map>
#include <string>

namespace rtengine {

constexpr double MINTEMP = 1500.0;
constexpr double MAXTEMP = 60000.0;
constexpr double MINGREEN = 0.02;
constexpr double MAXGREEN = 10.0;
constexpr double MINEQUAL = 0.8;
constexpr double MAXEQUAL = 1.5;
constexpr double INITIALBLACKBODY = 4000.0;


class ColorTemp {

private:
    double temp;
    double green;
    double equal;
    std::string method;
    static void clip(double &temp, double &green);
    static void clip(double &temp, double &green, double &equal);
    void temp2mul(double temp, double green, double equal, double& rmul, double& gmul, double& bmul) const;

public:
    ColorTemp() : temp(-1.), green(-1.), equal (1.), method("Custom") {}
    explicit ColorTemp(double e) : temp(-1.), green(-1.), equal (e), method("Custom") {}
    ColorTemp(double t, double g, double e, const std::string &m);
    ColorTemp(double mulr, double mulg, double mulb, double e);

    void update(const double rmul, const double gmul, const double bmul, const double equal)
    {
        this->equal = equal;
        mul2temp (rmul, gmul, bmul, this->equal, temp, green);
    }
    
    void useDefaults(const double equal)
    {
        temp = 6504;    // Values copied from procparams.cc
        green = 1.0;
        this->equal = equal;
    }

    inline std::string getMethod() const
    {
        return method;
    }
    
    inline double getTemp() const
    {
        return temp;
    }
    
    inline double getGreen() const
    {
        return green;
    }
    
    inline double getEqual() const
    {
        return equal;
    }

    void  getMultipliers(double &mulr, double &mulg, double &mulb) const
    {
        temp2mul(temp, green, equal, mulr, mulg, mulb);
    }

    void mul2temp(const double rmul, const double gmul, const double bmul, const double equal, double& temp, double& green) const;

    bool operator==(const ColorTemp& other) const
    {
        return fabs(temp - other.temp) < 1e-10 && fabs(green - other.green) < 1e-10;
    }

    bool operator!=(const ColorTemp& other) const
    {
        return !(*this == other);
    }


    static void temp2xy(double temp, double &x, double &y);
};

} // namespace rtengine
#endif
