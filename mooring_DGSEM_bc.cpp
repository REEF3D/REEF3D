/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"mooring_DGSEM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void mooring_DGSEM::updateBC(lexer* p)
{
	// Inlet Dirichlet bc
	
	r_xI = xs;
	r_yI = ys;
	r_zI = zs;
	v_xI = 0.0;
	v_yI = 0.0;
	v_zI = 0.0;	
	
	
	// Outlet Dirichlet bc
	
	// Coupling to 6DOF
/*	
	double xe = p->X311_xe[line];
	double ye = p->X311_ye[line];
	double ze = p->X311_ze[line];
*/
	
	// Ramp up
	double rampTime = 4.0;
	
	double Q = (t_mooring < rampTime) ? 0.5*(sin(1.0/rampTime*PI*t_mooring-PI/2)+1.0) : 1.0;	
	
	// Predefined path
/*
	double xe = 32.554 + Q*0.5*t_mooring;
	double ye = 0.0;
	double ze = 3.3 + Q*0.5*t_mooring;
*/

	// Circular motion
	double xe = 32.554;// + Q*0.2*cos(2*PI/3.5*t_mooring);
	double ye = 0.0;
	double ze = 3.3;// + Q*0.2*sin(2*PI/3.5*t_mooring);
	
	
	
	double timediff = p->simtime + p->dt - t_mooring_n;
	
	dtau += dt;
	
	a_xO = 2.0/timediff*(xe - r_xOn - v_xOn*timediff);
	a_yO = 2.0/timediff*(ye - r_yOn - v_yOn*timediff);
	a_zO = 2.0/timediff*(ze - r_zOn - v_zOn*timediff);

	r_xO = r_xOn + v_xOn*dtau + a_xO/2.0*dtau*dtau;
	r_yO = r_yOn + v_yOn*dtau + a_yO/2.0*dtau*dtau;
	r_zO = r_zOn + v_zOn*dtau + a_zO/2.0*dtau*dtau;
	
	v_xO = v_xOn + a_xO*dtau;
	v_yO = v_yOn + a_yO*dtau;
	v_zO = v_zOn + a_zO*dtau;	

	
/*
	// Inlet Neumann bc
	double t_xI, t_yI, t_zI, TI;
	t_xI = 1.0;
	t_yI = 0.0;
	t_zI = 0.0;
	TI = 0.0;
	
	neumann_bc
	(
		q_x_[0][0],q_y_[0][0],q_z_[0][0],
		v_x_[0][0],v_y_[0][0],v_z_[0][0],
		t_xI,t_yI,t_zI,TI,
		-1.0,0,0,0
	);
*/
/*
	// Outlet Neumann bc
	double t_xO, t_yO, t_zO, TO;
	t_xO = 1.0;
	t_yO = 0.0;
	t_zO = 0.0;
	TO = 2000.0;
	
	neumann_bc
	(
		q_x_[H-1][P],q_y_[H-1][P],q_z_[H-1][P],
		v_x_[H-1][P],v_y_[H-1][P],v_z_[H-1][P],
		t_xO,t_yO,t_zO,TO,
		1.0,H-1,1,P
	);
*/
}

void mooring_DGSEM::getBoundaryFluxes
(
	lexer *p, 
	double **& r_x_, double **& r_y_, double **& r_z_, 
	double **& q_x_, double **& q_y_, double **& q_z_, 
	double **& v_x_, double **& v_y_, double **& v_z_
)
{
	dirichlet_bc
	(
		r_x_[0][0],r_y_[0][0],r_z_[0][0],
		v_x_[0][0],v_y_[0][0],v_z_[0][0],
		r_xI,r_yI,r_zI,
		v_xI,v_yI,v_zI,
		-1.0,
		0,0,0,0
	);

	dirichlet_bc
	(
		r_x_[H-1][P],r_y_[H-1][P],r_z_[H-1][P],
		v_x_[H-1][P],v_y_[H-1][P],v_z_[H-1][P],
		r_xO,r_yO,r_zO,
		v_xO,v_yO,v_zO,
		1.0,
		H-1,H,1,P
	);		
}

void mooring_DGSEM::dirichlet_bc
(
	double& r_x_ij, double& r_y_ij, double& r_z_ij, 
	double& v_x_ij, double& v_y_ij, double& v_z_ij,
	double& r_xB, double& r_yB, double& r_zB, 
	double& v_xB, double& v_yB, double& v_zB, 
	double nx,
	int hi, int hf, int pf, int pj
)
{
	// Calculate boundary interface fluxes
	dfr_x[hi][pf] = - lambda_ct[hf]/2*(r_x_ij - r_xB);
	dfr_y[hi][pf] = - lambda_ct[hf]/2*(r_y_ij - r_yB);
	dfr_z[hi][pf] = - lambda_ct[hf]/2*(r_z_ij - r_zB);

	dfq_x[hi][pf] = nx/2.0*(fq_x[hi][pj] + v_xB);
	dfq_y[hi][pf] = nx/2.0*(fq_y[hi][pj] + v_xB);
	dfq_z[hi][pf] = nx/2.0*(fq_z[hi][pj] + v_xB);

	dfv_x[hi][pf] = - lambda_ct[hf]/2*(v_x_ij - v_xB*gamma);
	dfv_y[hi][pf] = - lambda_ct[hf]/2*(v_y_ij - v_yB*gamma);
	dfv_z[hi][pf] = - lambda_ct[hf]/2*(v_z_ij - v_zB*gamma);		
}


void mooring_DGSEM::neumann_bc
(
	double& q_x_ij, double& q_y_ij, double& q_z_ij, 
	double& v_x_ij, double& v_y_ij, double& v_z_ij,
	double& t_xB, double& t_yB, double& t_zB, double& TB,
	double nx, 
	int hi, int pf, int pj
)
{
	double cn = sqrt(T[hi][pj]/(qmag[hi][pj]*gamma));
	
	// Calculate boundary interface fluxes
	dfr_x[hi][pf] = 0.0;
	dfr_y[hi][pf] = 0.0;
	dfr_z[hi][pf] = 0.0;

	dfq_x[hi][pf] = nx*(fq_x[hi][pj] + v_x_ij/gamma)/2.0 - cn/2*(q_x_ij - qmag[hi][pj]*t_xB);
	dfq_y[hi][pf] = nx*(fq_y[hi][pj] + v_y_ij/gamma)/2.0 - cn/2*(q_y_ij - qmag[hi][pj]*t_yB);
	dfq_z[hi][pf] = nx*(fq_z[hi][pj] + v_z_ij/gamma)/2.0 - cn/2*(q_z_ij - qmag[hi][pj]*t_zB);

	dfv_x[hi][pf] = nx*(fv_x[hi][pj] + TB*t_xB)/2.0;
	dfv_y[hi][pf] = nx*(fv_y[hi][pj] + TB*t_yB)/2.0;
	dfv_z[hi][pf] = nx*(fv_z[hi][pj] + TB*t_zB)/2.0;	
}
