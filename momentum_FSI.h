/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
class sixdof_fsi;
class net;

using namespace std;

#ifndef MOMENTUM_FSI_H_
#define MOMENTUM_FSI_H_

class momentum_fsi : public momentum, public bcmom
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	momentum_fsi(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
	virtual ~momentum_fsi();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*);
	virtual void predictor(lexer*, fdm*, ghostcell*, momentum*, vrans*);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);
    virtual void fillaij1(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij2(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij3(lexer*, fdm*, ghostcell*, solver*);


	void starti(lexer*, fdm*, ghostcell*, sixdof_fsi*, vrans*, vector<net*>&);
    void ini(lexer*, fdm*, ghostcell*, sixdof_fsi*, vrans*, vector<net*>&);

    void updateFlags(lexer*, fdm*, ghostcell*);
    void predictor(lexer*, fdm*, ghostcell*, momentum*);
    void forcing(lexer*, fdm*, ghostcell*, sixdof_fsi*,field&,field&,field&,field&,field&,field&,double,vrans*,vector<net*>&);

    void fieldExt1(lexer*, fdm*, ghostcell*);
    void fieldExt2(lexer*, fdm*, ghostcell*);    
    void fieldExt3(lexer*, fdm*, ghostcell*);
   
    void forces(lexer*, fdm*, ghostcell*, sixdof_fsi*, field&,field&,field&,bool, double);
    void forces_surface(lexer*, fdm*, ghostcell*, sixdof_fsi*, bool);


    double Xfb, Yfb, Zfb, Mfb, Nfb, Kfb, cd, cq, cl;

private:

    void forcing1(lexer*, fdm*, ghostcell*);
    void forcing2(lexer*, fdm*, ghostcell*);
    void forcing3(lexer*, fdm*, ghostcell*);

    void forcing_uf(lexer*, fdm*, ghostcell*);
    void forcing_vf(lexer*, fdm*, ghostcell*);
    void forcing_wf(lexer*, fdm*, ghostcell*);
	
    void iniStencil1(lexer*, fdm*, ghostcell*);
    void iniStencil2(lexer*, fdm*, ghostcell*);
    void iniStencil3(lexer*, fdm*, ghostcell*);
    Eigen::VectorXd polBasis(const double&, const double&, const double&);
    double powPol(const double&, const double&);

    void updateStencil1(lexer*, fdm*, ghostcell*);
    void updateStencil2(lexer*, fdm*, ghostcell*);
    void updateStencil3(lexer*, fdm*, ghostcell*);
    
    void updateStencil_x(lexer*, fdm*, ghostcell*);
    void updateStencil_y(lexer*, fdm*, ghostcell*);
    void updateStencil_z(lexer*, fdm*, ghostcell*);
    void updateStencil_p(lexer*, fdm*, ghostcell*);
    
    void collectStencil
    (
        Eigen::Matrix4d&, Eigen::Matrix3d&,
        field&,
        const int, const int, const int,
        const double&, const double&, const double&,
        const double&, const double&, const double&,
        const double&, const double&, const double&,
        const double&, const double&, const double&,
        const double&, const double&, const double&
    );

    void findBoundPoint(const Eigen::MatrixXd&, double&, double&, double&);
    
    void reconstruct_pressure(lexer*, fdm*, ghostcell*);
    Eigen::MatrixXd updateRBDField_2D(lexer*, fdm*, ghostcell*);
    Eigen::MatrixXd updateRBDField_3D(lexer*, fdm*, ghostcell*);

    double Hsolidface(lexer*, fdm*, int, int, int);
    
    field1 urk1, urk2, un, uf, fx, flagx, gradPx;
	field2 vrk1, vrk2, vn, vf, fy, flagy, gradPy;
	field3 wrk1, wrk2, wn, wf, fz, flagz, gradPz;
    field4 flagp;

	convection *pconvec;
	diffusion *pdiff, *pdiff_e;
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
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);    
    
	double starttime;
};

#endif

