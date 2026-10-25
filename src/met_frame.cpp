/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"met_frame.h"
#include<cmath>

namespace
{
    const double pi = 3.14159265358979323846;
    const double rad = pi/180.0;
}

met_frame::met_frame(double lon0_, double lat0_, double rot_, double X0_, double Y0_)
{
    set = true;
    lon0 = lon0_;
    lat0 = lat0_;
    rot  = rot_;
    X0   = X0_;
    Y0   = Y0_;

    lam0 = lon0*rad;
    phi0 = lat0*rad;
    c = std::cos(rot*rad);
    s = std::sin(rot*rad);
}

void met_frame::to_model(double lon, double lat, double &x, double &y) const
{
    // spherical transverse Mercator about lam0, northing from phi0
    const double lam = lon*rad - lam0;
    const double phi = lat*rad;
    const double b = std::cos(phi)*std::sin(lam);
    const double e = 0.5*R*std::log((1.0 + b)/(1.0 - b));
    const double n = R*(std::atan2(std::tan(phi), std::cos(lam)) - phi0);

    x = c*e + s*n + X0;
    y = -s*e + c*n + Y0;
}

void met_frame::to_geo(double x, double y, double &lon, double &lat) const
{
    const double xx = x - X0;
    const double yy = y - Y0;
    const double e = c*xx - s*yy;
    const double n = s*xx + c*yy;
    const double d = n/R + phi0;

    lat = std::asin(std::sin(d)/std::cosh(e/R))/rad;
    lon = (lam0 + std::atan2(std::sinh(e/R), std::cos(d)))/rad;
}

double met_frame::north(double x, double y) const
{
    double lon, lat, x1, y1;

    to_geo(x,y,lon,lat);
    to_model(lon,lat+1.0e-4,x1,y1);

    double x0, y0;
    to_model(lon,lat,x0,y0);

    return std::atan2(y1-y0, x1-x0);
}

double met_frame::north_geo(double lon, double lat) const
{
    double x0, y0, x1, y1;

    to_model(lon,lat,x0,y0);
    to_model(lon,lat+1.0e-4,x1,y1);

    return std::atan2(y1-y0, x1-x0);
}
