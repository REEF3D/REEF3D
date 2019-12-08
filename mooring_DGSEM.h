/*--------------------------------------------------------------------
REEF3D
Copyright 2019 Tobias Martin

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"mooring.h"
#include"mooring_Catenary.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field5.h"
#include"fieldint5.h"
#include"vec.h"
#include<fstream>
#include<iostream>
#include<vector>


using namespace std;

#ifndef MOORING_DGSEM_H_
#define MOORING_DGSEM_H_

class mooring_DGSEM : public mooring
{
public:
	mooring_DGSEM(int);
	virtual ~mooring_DGSEM();
	
	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void initialize(lexer*, fdm*, ghostcell*);
	virtual void mooringForces(double&, double&, double&);
	
private:

	// Initialisation
	void ini_parallel(lexer*, fdm*, ghostcell*);
	void facenodegen(lexer*, double, double);
	void getqp(lexer*);
	void getDr(lexer*);
	void jacobiP(double, double, int, double*&);
	void getV(lexer*);
	void getsurfInt(lexer*);
	void nodegen(lexer*);
	void iniCond(lexer*, fdm*, ghostcell*);

	// Runtime
	void updateFluidVel(lexer*, fdm*, ghostcell*, int);
	void updateFields(lexer*, fdm*, ghostcell*);
	void updateBC(lexer*);
	void saveMooringPoint(lexer*);
	void startRK3TVD(lexer*, fdm*, ghostcell*);
	void getL(lexer*, double**&, double**&, double**&, double**&, double**&, double**&, double**&, double**&, double**&);
	double Drf(double*&, int);
	double sIntdf(double*&, int);
	void getBoundaryFluxes(lexer*, double**&, double**&, double**&, double**&, double**&, double**&, double**&, double**&, double**&);
	void dirichlet_bc(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double, int, int, int, int);
	void neumann_bc(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double, int, int, int);
	void getForces(double**&, double**&, double**&, double**&);	
	void limitFields(lexer*, double**&, double**&, double**&, double**&, double**&, double**&, double**&, double**&, double**&);
	void startLimit(lexer*, double**&);
	double minMod(const double&, const double&, const double&);
	void slopeLimit(lexer *, int, double*&, double*&, double*&);

	// Print
	void print(lexer*);
	
	// ------ 
	
	// Initialisation
	mooring_Catenary *pcatenary;
	
	// Parallelisation
	int line;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
	// Material constants
	double gamma, EA, L, rho_c, d_c;
	
	// Mesh
	int P, H;
	double rx;
	double *vx, *r;
	double **Dr, **V, **invV, **x, **sInt;
	double xs, ys, zs;

	// Fields and fluxes
	double **T, **qmag, **t_x, **t_y, **t_z;
	double **Fx, **Fy, **Fz;
	double **r_x, **r_y, **r_z, **q_x, **q_y, **q_z, **v_x, **v_y, **v_z, **v_xn, **v_yn, **v_zn;
	double **r_x1, **r_y1, **r_z1, **q_x1, **q_y1, **q_z1, **v_x1, **v_y1, **v_z1;
	double **r_x2, **r_y2, **r_z2, **q_x2, **q_y2, **q_z2, **v_x2, **v_y2, **v_z2;
	double **Lr_x, **Lr_y, **Lr_z, **Lq_x, **Lq_y, **Lq_z, **Lv_x, **Lv_y, **Lv_z;
	double **fq_x, **fq_y, **fq_z, **fv_x, **fv_y, **fv_z;
	double **dfr_x, **dfr_y, **dfr_z, **dfq_x, **dfq_y, **dfq_z, **dfv_x, **dfv_y, **dfv_z;
	double *lambda_cn, *lambda_ct;
	double ***vf, ***vfn, ***af;	
		
	// Boundary conditions
	double r_xI, r_yI, r_zI, v_xI, v_yI, v_zI;
	double r_xO, r_yO, r_zO, v_xO, v_yO, v_zO;
	double r_xOn, r_yOn, r_zOn, v_xOn, v_yOn, v_zOn, a_xO, a_yO, a_zO;
	
	// Time control
	double cfl, relFac, dt, dtm, dtau, t_mooring_n, t_mooring;
	
	// Forces
	double Xme_, Yme_, Zme_;
	
	// Print
	char name[100];
	ofstream eTout;
	double printtime;
	
};

#endif
