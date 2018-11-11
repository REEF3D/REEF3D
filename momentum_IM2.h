/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"momentum.h"
#include"ibcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"

class discrete;
class diffusion;
class pressure;
class turbulence;
class solver;
class poisson;

using namespace std;

#ifndef MOMENTUM_IM2_H_
#define MOMENTUM_IM2_H_

class momentum_IM2 : public momentum, public ibcmom
{
public:
	momentum_IM2(lexer*, fdm*, ghostcell*, discrete*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
	virtual ~momentum_IM2();
	virtual void start(lexer*, fdm*, ghostcell*, momentum*);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);
    virtual void fillaij1(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij2(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij3(lexer*, fdm*, ghostcell*, solver*);

	void usource(lexer*,fdm*);
	void vsource(lexer*,fdm*);
	void wsource(lexer*,fdm*);
	field1 un,unn;
	field2 vn,vnn;
	field3 wn,wnn;

private:
	void irhs(lexer*,fdm*);
	void jrhs(lexer*,fdm*);
	void krhs(lexer*,fdm*);
	
    void clearrhs(lexer*,fdm*);
	int gcval_u, gcval_v, gcval_w;
	int count,q;
	double starttime;

	discrete *pdisc;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
};

#endif


