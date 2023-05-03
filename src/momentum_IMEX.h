/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"momentum.h"
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

using namespace std;

#ifndef MOMENTUM_IMEX_H_
#define MOMENTUM_IMEX_H_

class momentum_IMEX : public momentum, public bcmom
{
public:
	momentum_IMEX(lexer*, fdm*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
	virtual ~momentum_IMEX();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);

private:
    
    void eval_ex_F(lexer*,fdm*,ghostcell*,vrans*,field&,field&,field&,field&);
	void eval_ex_G(lexer*,fdm*,ghostcell*,vrans*,field&,field&,field&,field&);
	void eval_ex_H(lexer*,fdm*,ghostcell*,vrans*,field&,field&,field&,field&);

    field1 un, F0_ex, F1_ex, F2_ex, F1_im, F2_im;
    field2 vn, G0_ex, G1_ex, G2_ex, G1_im, G2_im;
    field3 wn, H0_ex, H1_ex, H2_ex, H1_im, H2_im;
    
    int gcval_u, gcval_v, gcval_w;
	double starttime, gamma, a11,a21,a22,ahat10,ahat20,ahat21,b1,b2,bhat1,bhat2,twoD;

	convection *pconvec;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
    fluid_update *pupdate;
};

#endif
