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

#ifndef MET_PROVIDER_H_
#define MET_PROVIDER_H_

#include"met_frame.h"
#include<string>

class lexer;
class ghostcell;
class met_source;

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200): the met provider of the rank.

One provider per process, shared by all grids of the rank (level 0, the
refined patches, the blocks): it owns the source (W 200 1 uniform,
W 200 2 met file), the model start time and the frame, and advances the
source once per model time. Every grid has its own met_forcing with the
fields of its cells; it registers the box of its cells (window) and
evaluates the source per cell (at).

  W 203 YYYYMMDD.HHMMSS  model time 0; 0: the first record of the file
  W 204 T                ramp of the forcing over T s, 0.5 (1 - cos(pi t/T))
  W 205 lon0 lat0 rot    frame, W 206 X0 Y0 (met_frame)
  W 207 law              drag law (met_drag.h)

The provider holds no field of a grid, only the source (static data
of the rank, Regression/unit/static_allowlist.txt).
--------------------------------------------------------------------*/

class met_provider
{
public:
    static met_provider* get(lexer*, ghostcell*);

    void window(double xs, double xe, double ys, double ye);
    void update(lexer*, ghostcell*, double simtime);
    void at(double x, double y, double &u, double &v, double &pa) const;
    double ramp(lexer*, double simtime) const;

    // checks of the W 200 input, from the driver for every model
    static void check(lexer*, ghostcell*);

    double tref = 0.0;      // absolute time of model time 0 [s since 1970]
    met_frame frame;

private:
    met_provider(lexer*, ghostcell*);
    void stop(lexer*, ghostcell*, const std::string&) const;

    static met_provider *instance;

    met_source *src = nullptr;
    bool has_t = false;
    double tnow = 0.0;
};

#endif
