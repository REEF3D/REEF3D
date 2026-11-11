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

#ifndef SFLOW_MOMENTUM_H_
#define SFLOW_MOMENTUM_H_

class lexer;
class fdm2D;
class ghostcell;
class sflow_pressure;
class solver2D;
class ioflow;
class sixdof;

#include<functional>

using namespace std;

class sflow_momentum 
{
public:
    virtual ~sflow_momentum() = default;

    // checkpoints (P 48, I 45): the arrays of the object that carry over to the next step
    virtual void ckpt_fields(const std::function<void(const char*, char, double*)>&, const std::function<void(const char*, char, int*)>&) {}
    // after a restart from a checkpoint (the state of a step: the start-up of the first step is done)
    virtual void ckpt_restored() {}
    // the water column was changed after the step (REEF3D::HYDRO losses): the semi-implicit scheme keeps its
    // own copy of the column for the next step (sflow_momentum_SI::si_restore)
    virtual void water_column(lexer*, fdm2D*) {}
    // the discharges over the open boundaries of the step [m3/s], inflow and outflow positive (REEF3D::HYDRO budget);
    // false: not known exactly (explicit schemes: the budget takes the ioflow discharges)
    virtual bool boundary_flow(double&, double&) { return false; }

	virtual void start(lexer*, fdm2D*, ghostcell*)=0;
    virtual void ini(lexer*, fdm2D*, ghostcell*)=0;
};

#endif
