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

#ifndef SFLOW_RHEOLOGY_H_
#define SFLOW_RHEOLOGY_H_

class lexer;
class fdm2D;
class slice;
class ghostcell;

using namespace std;

class sflow_rheology
{

public:

    virtual ~sflow_rheology() = default;

    // bed shear stress for the conserved momenta: args U, V, WL (cell centres)
	virtual void u_source(lexer*, fdm2D*, slice&, slice&, slice&)=0;
    virtual void v_source(lexer*, fdm2D*, slice&, slice&, slice&)=0;
    
    // rheology core (W 120): initial state, and the implicit basal friction step on the stage
    // momenta UH, VH with the stage depth WL (args: UH, VH, WL, alpha of the stage, last stage)
    virtual void ini(lexer*, fdm2D*, ghostcell*) {}
    virtual void implicit(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double, bool) {}

    // semi-implicit SFLOW (A 210 11/12/13): the basal friction as a factor on the momentum of the step,
    // s/Ustar of the backward-Euler friction step for the frictionless estimate Ustar of the new speed
    // (args: cell i, j; dir 0 the cell, 1 its face i+1/2, 2 its face j+1/2; depth, Ustar, dt), and the end
    // of the step (pore pressure; args UH, VH, WL of t^{n+1})
    virtual bool si_friction() const { return false; }
    virtual double si_factor(lexer*, fdm2D*, int, int, int, double, double, double) { return 1.0; }
    virtual void si_end(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&) {}

};

#endif
