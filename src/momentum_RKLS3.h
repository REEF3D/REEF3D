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
#include"momentum_forcing.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
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
class sixdof_base;
class net;
class fsi;

using namespace std;

#ifndef MOMENTUM_RKLS3_H_
#define MOMENTUM_RKLS3_H_

class momentum_RKLS3 : public momentum, public momentum_forcing, public bcmom
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	momentum_RKLS3(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*, fsi*);
	virtual ~momentum_RKLS3();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof_df_base*,vector<net*>&);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);
    virtual void fillaij1(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij2(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij3(lexer*, fdm*, ghostcell*, solver*);


private:

    void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);    
    
    field1 urk, Cu, Du, fx;
	field2 vrk, Cv, Dv, fy;
	field3 wrk, Cw, Dw, fz;

	convection *pconvec;
	diffusion *pdiff;
	diffusion *pdiff_e;
	pressure *ppress;
	poisson *ppois;
	density *pdensity;
    turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow; 
    fsi *pfsi;   
    
	int gcval_u, gcval_v, gcval_w;

    Eigen::Vector3d alpha, gamma, zeta;

	double starttime;
};

#endif

