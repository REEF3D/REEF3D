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

#ifndef SEDIMENT_INFLOW_H_
#define SEDIMENT_INFLOW_H_

#include<vector>

class lexer;

/*--------------------------------------------------------------------
Sediment supply at the inflow (S 129), SFLOW, NHFLOW and CFD
  S 129 0  equilibrium: bed load zero gradient at the inflow (the import
           of the first cell is its own transport), concentration S 124
  S 129 1  constant: bed load S 130 [m^2/s, sediment volume per width],
           concentration S 124
  S 129 2  time series sediment-inflow.dat: lines  t [s]  qb [m^2/s]  C [-],
           linear in time, the first/last value before/after the series;
           qb < 0: zero gradient, C < 0: equilibrium
The inflow concentration (S 124 or the time series), volume concentration:
  SFLOW S 120: the depth-averaged concentration of the inflow faces,
               < 0 the equilibrium concentration C_eq
  NHFLOW S 12: the depth-averaged concentration of the inflow ghost
               columns, distributed over the depth with the profile of
               the first column; < 0 zero gradient (the profile of the
               first column)
The bed load enters the Exner step as the transport of the inflow ghost
cells (the import of the first cell, in the direction of the flow); with
a graded bed (S 51) it is split with the fractions F_a,k of S 51.
--------------------------------------------------------------------*/

class sediment_inflow
{
public:
    sediment_inflow(lexer*);
    
    // the supply at the time p->simtime; < 0: equilibrium (zero gradient)
    double qb(lexer*) const;
    double conc(lexer*) const;
    
    bool bedload_on() const {return mode>0;}
    
private:
    double ipol(lexer*, int col) const;
    
    int mode;
    std::vector<double> t,q,c;
};

#endif
