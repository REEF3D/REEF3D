/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#ifndef MOMENTUM_FSI_H_
#define MOMENTUM_FSI_H_

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

    void forcing(lexer*, fdm*, ghostcell*, sixdof_df*,field&,field&,field&,field&,field&,field&,double,vrans*,vector<net*>&);
    double Hsolidface(lexer*, fdm*, int, int, int);
	
    void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);    
    
    field1 urk1, urk2, un, uf, fx, flagx, gradPx;
	field2 vrk1, vrk2, vn, vf, fy, flagy, gradPy;
	field3 wrk1, wrk2, wn, wf, fz, flagz, gradPz;
    field4 flagp;

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

    int r, r_y, radius, radius_y, K, colSize, rowSize;

    std::vector<Eigen::Matrix4d> stencil_x, stencil_y, stencil_z, stencil_p;
    std::vector<Eigen::Matrix3d> stencil_x_b, stencil_y_b, stencil_z_b, stencil_p_b;
    std::vector<Eigen::Array4d> surf_press;
    
    std::vector<Eigen::MatrixXd> stencil1, stencil2, stencil3;
    
	double starttime;
};

#endif

