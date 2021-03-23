/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"momentum.h"
#include"bcmom.h"
#include"diffusion.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include<vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class density;
class poisson;
class sixdof_df;
class net;

using namespace std;

#ifndef MOMENTUM_RK3_DF_H_
#define MOMENTUM_RK3_DF_H_

class momentum_RK3_df : public momentum, public bcmom
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	momentum_RK3_df(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
	virtual ~momentum_RK3_df();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);
    virtual void fillaij1(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij2(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij3(lexer*, fdm*, ghostcell*, solver*);

	void starti(lexer*, fdm*, ghostcell*, sixdof_df*, vrans*, vector<net*>&);

private:

    double Hsolidface(lexer*, fdm*, int, int, int);
	
    void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);    
    
    field1 udiff, urk1, urk2, fx;
	field2 vdiff, vrk1, vrk2, fy;
	field3 wdiff, wrk1, wrk2, fz;

	convection *pconvec;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	density *pdensity;
    turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;    
    
	int gcval_u, gcval_v, gcval_w, gcval_urk, gcval_vrk, gcval_wrk;

	double starttime;
};

#endif

