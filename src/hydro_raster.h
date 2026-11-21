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

#ifndef HYDRO_RASTER_H_
#define HYDRO_RASTER_H_

#include<string>
#include<vector>

/*--------------------------------------------------------------------
REEF3D::HYDRO: reader of an ESRI ASCII grid (.asc), the raster format of
the HYDRO inputs (rain R 10 3, Manning n R 50). No lexer or MPI
dependency.

  ncols         nc
  nrows         nr
  xllcorner     x0        (or xllcenter: the centre of the lower left cell)
  yllcorner     y0        (or yllcenter)
  cellsize      cs        (or dx and dy)
  NODATA_value  nd        (optional)
  values, nr rows of nc values, the first row is the northern one (largest y)

Coordinates are world coordinates (as geo.dat, before the rotation and
offset of the grid, lexer::XYin). A point takes the value of the raster
cell it lies in (no interpolation: rain totals and land-use classes are
cell values); outside the raster and at NODATA, at() returns false.
--------------------------------------------------------------------*/

class hydro_raster
{
public:
    bool read(const std::string &file, std::string &error);
    bool at(double x, double y, double &val) const;

    int nc = 0, nr = 0;
    double x0 = 0.0, y0 = 0.0, dx = 1.0, dy = 1.0, nodata = -9999.0;
    bool has_nodata = false;

private:
    std::vector<double> v;      // row-major from the south: v[jr*nc + ir]
};

#endif
