/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sflow_rheology.h"
#include"increment.h"
#include"slice4.h"

#ifndef SFLOW_RHEOLOGY_F_H_
#define SFLOW_RHEOLOGY_F_H_

#define HXIJ (fabs(b->hx(i,j))>1.0e-10?b->hx(i,j):1.0e20)
#define HYIJ (fabs(b->hy(i,j))>1.0e-10?b->hy(i,j):1.0e20)
#define HPIJ (fabs(b->hp(i,j))>1.0e-10?b->hp(i,j):1.0e20)

using namespace std;

class sflow_rheology_f : public sflow_rheology, public increment
{

public:
    sflow_rheology_f(lexer*);
	virtual ~sflow_rheology_f();
    
	virtual void u_source(lexer*, fdm2D*, slice&, slice&);
    virtual void v_source(lexer*, fdm2D*, slice&, slice&);

private:
    double bingham(lexer*, fdm2D*, double, double, double, double);
    
    double cf,manning;
    
    double tau_zx,tau_zy,u_abs;
    double tau0,val;
    double press;
    double tanphi;
    double hc,denom;
};

#endif
