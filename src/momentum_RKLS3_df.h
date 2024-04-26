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
class sixdof;
class net;
class fsi;

using namespace std;

#ifndef MOMENTUM_RKLS3_DF_H_
#define MOMENTUM_RKLS3_DF_H_

class momentum_RKLS3_df : public momentum, public bcmom
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	momentum_RKLS3_df(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
	virtual ~momentum_RKLS3_df();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof*,vector<net*>&);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);
    virtual void fillaij1(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij2(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij3(lexer*, fdm*, ghostcell*, solver*);

	void starti(lexer*, fdm*, ghostcell*, sixdof*, vrans*, vector<net*>&, fsi*);

private:

    double Hsolidface(lexer*, fdm*, int, int, int);
	
    void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);    
    
    field1 urk, Cu, Du, fx;
	field2 vrk, Cv, Dv, fy;
	field3 wrk, Cw, Dw, fz;

	convection *pconvec;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
    turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;    
    
	int gcval_u, gcval_v, gcval_w;

    Eigen::Vector3d alpha, gamma, zeta;

	double starttime;
};

#endif

