/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef MOMENTUM_RK3_V2_H_
#define MOMENTUM_RK3_V2_H_

#include"momentum.h"
#include"momentum_forcing.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class poisson;
class fluid_update;
class nhflow;
class sixdof;
class fsi;

using namespace std;

class momentum_RK3_v2 : public momentum, public momentum_forcing, public bcmom
{
public:
	momentum_RK3_v2(lexer*, fdm*, convection*, diffusion*, pressure*, poisson*, 
                turbulence*, solver*, solver*, ioflow*, fsi*);
	virtual ~momentum_RK3_v2();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof*);

    field1 urk1,urk2,Frk1,Frk2,fx;
	field2 vrk1,vrk2,Grk1,Grk2,fy;
	field3 wrk1,wrk2,Hrk1,Hrk2,fz;

private:
    fluid_update *pupdate;
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
    
	int gcval_u, gcval_v, gcval_w;
	double starttime;

	convection *pconvec;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
    nhflow *pnh;
    fsi *pfsi;
};

#endif
