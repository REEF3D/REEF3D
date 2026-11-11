/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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


#include"nhflow_gradient.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"nhflow_curv.h"
#include<cmath>

nhflow_gradient::nhflow_gradient(lexer* pp) : tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),dx(pp->DXM)
{
    p=pp;
    
    grad=0.0;
}

nhflow_gradient::~nhflow_gradient()
{

}

// **********************************************************
// wall-function bed (A 518 1)
// **********************************************************
// The ghost cells of the wall-function bed are mirrored (zero gradient, gc_startV.cpp): the bed shear is the
// wall friction of A 519 (nhflow_bcmom). The vertical derivative of the horizontal velocity in the bed cell
// is then the one of the log law the friction is built on (nhflow_bed_drag):
//    du/dz = u_tau/(kappa z) along u,  u_tau = |u|/u+,  u+ = (1/kappa) max(1, ln(30 z/ks)),
//    z = 0.5 dz of the cell (>= ks/30)   ->   du/dz = u/(kappa z u+)   (component-wise, linear in u)
// It enters the strain rate of the bed cell: the eddy viscosity of LES, the limiters of the RANS models
// (their production in the bed cell is the wall function's, nhflow_strain::wallf_update) and the vorticity.
// The ghost value 0 of A 518 2 gave du/dz = u/dz, ln(30 z/ks)/2 times the log law (3-4 times for
// z/ks = 10-100), the mirrored ghosts alone give 0.
// Only for U and V of the fdm registered with bedwall_set (the strain of the turbulence models): the
// derivatives are used for other fields too (stresses, W, scalars), these keep the ghost cells.

void nhflow_gradient::bedwall_set(fdm_nhf *d)
{
    bw_d = d;
}

bool nhflow_gradient::bedwall(const double *f)
{
    return bw_d!=nullptr && k==0 && p->A518==1 && p->A519>0 && p->flag4(i,j,k-1)<0
        && (f==bw_d->U || f==bw_d->V);
}

// d f/d sigma of the bed cell from the log law (bedwall() true)
double nhflow_gradient::dfdsig_bed(const double *f)
{
    const double kappa=0.4;
    const double WLval = bw_d->WL(i,j);
    double ks = bw_d->ks(i,j);
    
    if(ks<=0.0)
    ks=0.0001;
    
    double z = 0.5*p->DZN[KP]*WLval;
    
    if(30.0*z<ks)
    z=ks/30.0;
    
    const double uplus = (1.0/kappa)*MAX(1.0,log(30.0*(z/ks)));
    
    // du/dz = u/(kappa z u+), d/dsigma = WL d/dz
    return WLval*f[IJK]/(kappa*z*uplus);
}

// d f/d sigma at the cell centre, central (the bed cell of the wall-function bed: log law)
double nhflow_gradient::dfdsig_c(const double *f)
{
    if(bedwall(f))
    return dfdsig_bed(f);
    
    return (f[IJKp1] - f[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
}

// curvilinear grid (G 50): the horizontal derivatives along the sigma layer by the chain rule
static inline double curv_dx(lexer *p, double *f, int i, int j, int k, int comp)
{
    double fx,fy;
    nhflow_curv::grad5(p,i,j,f[Ip1JK],f[Im1JK],f[IJp1K],f[IJm1K],fx,fy);
    return comp==1 ? fx : fy;
}

// **********************************************************
// DUXDX2
// **********************************************************

double nhflow_gradient::dudx(double *U)
{
	grad = (p->geo_curv==1 ? curv_dx(p,U,i,j,k,1) : (U[Ip1JK] - U[Im1JK])/(p->DXP[IP]+p->DXP[IM1])) 
    
         + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dfdsig_c(U);

	return grad;
}

double nhflow_gradient::dudy(double *U)
{
	grad = (p->geo_curv==1 ? curv_dx(p,U,i,j,k,2) : (U[IJp1K] - U[IJm1K])/(p->DYP[JP]+p->DYP[JM1]))
    
         + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dfdsig_c(U);

	return grad;
}

double nhflow_gradient::dudz(double *U)
{
	grad = p->sigz[IJ]*(U[IJKp1] - U[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
    
    if(k==0)
    grad = p->sigz[IJ]*(U[IJK] - U[IJKm1])/(p->DZP[KM1]);
    
    if(bedwall(U))
    grad = p->sigz[IJ]*dfdsig_bed(U);
    
    //if(k==p->knoz-1)
    //grad = p->sigz[IJ]*(U[IJK] - U[IJKm1])/(p->DZP[KM1]);

	return grad;
}

// **********************************************************
// DVDX2
// **********************************************************

double nhflow_gradient::dvdx(double *V)
{
	grad = (p->geo_curv==1 ? curv_dx(p,V,i,j,k,1) : (V[Ip1JK] - V[Im1JK])/(p->DXP[IP]+p->DXP[IM1]))
    
         + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dfdsig_c(V);

	return grad;
}

double nhflow_gradient::dvdy(double *V)
{
	grad = (p->geo_curv==1 ? curv_dx(p,V,i,j,k,2) : (V[IJp1K] - V[IJm1K])/(p->DYP[JP]+p->DYP[JM1]))
    
         + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dfdsig_c(V);

	return grad;
}

double nhflow_gradient::dvdz(double *V)
{
	grad = p->sigz[IJ]*(V[IJKp1] - V[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
    
    if(k==0)
    grad = p->sigz[IJ]*(V[IJK] - V[IJKm1])/(p->DZP[KM1]);
    
    if(bedwall(V))
    grad = p->sigz[IJ]*dfdsig_bed(V);
    
    //if(k==p->knoz-1)
    //grad = p->sigz[IJ]*(V[IJK] - V[IJKm1])/(p->DZP[KM1]);

	return grad;
}

// **********************************************************
// DZX2
// **********************************************************

double nhflow_gradient::dwdx(double *W)
{
	grad = (p->geo_curv==1 ? curv_dx(p,W,i,j,k,1) : (W[Ip1JK] - W[Im1JK])/(p->DXP[IP]+p->DXP[IM1]))
        
         + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(W[IJKp1] - W[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

double nhflow_gradient::dwdy(double *W)
{
	grad = (p->geo_curv==1 ? curv_dx(p,W,i,j,k,2) : (W[IJp1K] - W[IJm1K])/(p->DYP[JP]+p->DYP[JM1]))
        
         + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(W[IJKp1] - W[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}


double nhflow_gradient::dwdz(double *W)
{
	grad = p->sigz[IJ]*(W[IJKp1] - W[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
    
    if(k==0)
    grad = p->sigz[IJ]*(W[IJK] - W[IJKm1])/(p->DZP[KM1]);
    
    //if(k==p->knoz-1)
    //sgrad = p->sigz[IJ]*(W[IJK] - W[IJKm1])/(p->DZP[KM1]);

	return grad;
}


