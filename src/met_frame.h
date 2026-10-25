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

#ifndef MET_FRAME_H_
#define MET_FRAME_H_

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200): geographic frame of the model grid.
No lexer, unit test met_test.cpp.

W 205 lon0 lat0 rot   transverse Mercator on a sphere (R = 6371 km) about
                      lon0, lat0 [deg]; the model x axis points east,
                      rotated counter-clockwise by rot [deg]
W 206 X0 Y0           model coordinates of (lon0, lat0) [m], default 0 0

The same frame as tools/seastate_forcing.py --origin LON0 LAT0
--rotation ROT (X0 = Y0 = 0) and tools/met_forcing.py: the met files are
written on a grid in model coordinates, and the frame turns lon/lat
(storm tracks, latitude for the Coriolis parameter) into model
coordinates and back.

  to_model(lon,lat,x,y)   [deg] -> [m]
  to_geo(x,y,lon,lat)     [m] -> [deg]
  north(x,y)              direction of true north in the model frame
                          [rad, counter-clockwise from the model x axis]
  north_geo(lon,lat)      the same at a geographic point
--------------------------------------------------------------------*/

class met_frame
{
public:
    met_frame() = default;
    met_frame(double lon0, double lat0, double rot, double X0, double Y0);

    void to_model(double lon, double lat, double &x, double &y) const;
    void to_geo(double x, double y, double &lon, double &lat) const;
    double north(double x, double y) const;
    double north_geo(double lon, double lat) const;     // the same at (lon, lat)

    bool set = false;
    double lon0 = 0.0, lat0 = 0.0, rot = 0.0, X0 = 0.0, Y0 = 0.0;   // as given [deg], [m]

    static constexpr double R = 6371000.0;

private:
    double lam0 = 0.0, phi0 = 0.0, c = 1.0, s = 0.0;
};

#endif
