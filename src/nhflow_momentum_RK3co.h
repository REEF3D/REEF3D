/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"momentum.h"
#include"bcmom.h"
#include"field4.h"
#include"slice4.h"

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class poisson;
class fluid_update;
class nhflow;
class nhflow_fsf;

using namespace std;

#ifndef NHFLOW_MOMENTUM_RK3co_H_
#define NHFLOW_MOMENTUM_RK3co_H_

class nhflow_momentum_RK3co : public momentum, public bcmom
{
public:
	nhflow_momentum_RK3co(lexer*, fdm*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*, nhflow*, nhflow_fsf*);
	virtual ~nhflow_momentum_RK3co();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*);
    virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);

    field4 udiff,urk1,urk2;
	field4 vdiff,vrk1,vrk2;
	field4 wdiff,wrk1,wrk2;
    
    slice4 etark1,etark2;

private:
    fluid_update *pupdate;
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	
    void timecheck(lexer*,fdm*,ghostcell*,field&,field&,field&);
    
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
    nhflow_fsf *pnhfsf;
};

#endif
