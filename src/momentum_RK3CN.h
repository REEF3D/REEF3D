/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

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
class onephase;
class solver;
class poisson;
class fluid_update;
class nhflow;
class sixdof;
class fsi;

using namespace std;

#ifndef MOMENTUM_RK3CN_H_
#define MOMENTUM_RK3CN_H_

class momentum_RK3CN : public momentum, public momentum_forcing, public bcmom
{
public:
	momentum_RK3CN(lexer*, fdm*, convection*, diffusion*, pressure*, poisson*, turbulence*, onephase*, solver*, solver*, ioflow*, fsi*);
	virtual ~momentum_RK3CN();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof*,vector<net*>&);
        virtual void utimesave(lexer*, fdm*, ghostcell*);
        virtual void vtimesave(lexer*, fdm*, ghostcell*);
        virtual void wtimesave(lexer*, fdm*, ghostcell*);

    field1 udiff,urk1,urk2,fx;
	field2 vdiff,vrk1,vrk2,fy;
	field3 wdiff,wrk1,wrk2,fz;

private:
        fluid_update *pupdate;
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
        void addirhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
        void addjrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
        void addkrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	
        void timecheck(lexer*,fdm*,ghostcell*,field&,field&,field&);
    
	int gcval_u, gcval_v, gcval_w;
	double starttime;

	convection *pconvec;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
    onephase *poneph;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
    nhflow *pnh;
    fsi *pfsi;
};

#endif
