/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"convection.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"turbulence.h"
#include"solver.h"

using namespace std;

#ifndef MOMENTUM_AB2_H_
#define MOMENTUM_AB2_H_

class momentum_AB2 : public momentum, public bcmom
{
public:
	momentum_AB2(lexer*, fdm*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
	virtual ~momentum_AB2();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof_df_base*,vector<net*>&);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);

	field1 uab;
	field2 vab;
	field3 wab;

private:
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
};

#endif
