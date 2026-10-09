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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SFLOW_PRESSURE_NH_H_
#define SFLOW_PRESSURE_NH_H_

#include"sflow_pressure.h"

class sflow_amr;

using namespace std;

// non-hydrostatic pressure (A 220 1-3) split into its steps, so that the mesh refinement
// (sflow_amr) can solve one pressure on level 0 and the patches together:
//   assemble: rows and right-hand side of the grid (sflow_pjm_quad: the bed acceleration first)
//   correct:  velocity correction with the solved pressure
// the level-0 instance hands the solve to sflow_amr::nh_solve when patches exist

class sflow_pressure_nh : public sflow_pressure
{
public:
    sflow_amr *amr = nullptr;

    virtual void assemble(lexer*, fdm2D*, ghostcell*, slice &UH, slice &VH, slice &WL, slice &Un, slice &Vn, double alpha)=0;
    virtual void correct(lexer*, fdm2D*, slice &UH, slice &VH, slice &WH, slice &WL, double alpha)=0;
    virtual int is_active(lexer*, fdm2D*, int, int)=0;
    virtual int gcval() const=0;
    // coefficients of the bed pressure in the u, v correction and of the w correction (G 7 1:
    // the synchronisation projection of sflow_amr)
    virtual void coef(double &cbq, double &cwq) const=0;
};

#endif
