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
    
         + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(U[IJKp1] - U[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

double nhflow_gradient::dudy(double *U)
{
	grad = (p->geo_curv==1 ? curv_dx(p,U,i,j,k,2) : (U[IJp1K] - U[IJm1K])/(p->DYP[JP]+p->DYP[JM1]))
    
         + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(U[IJKp1] - U[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

double nhflow_gradient::dudz(double *U)
{
	grad = p->sigz[IJ]*(U[IJKp1] - U[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
    
    if(k==0)
    grad = p->sigz[IJ]*(U[IJK] - U[IJKm1])/(p->DZP[KM1]);
    
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
    
         + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(V[IJKp1] - V[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

double nhflow_gradient::dvdy(double *V)
{
	grad = (p->geo_curv==1 ? curv_dx(p,V,i,j,k,2) : (V[IJp1K] - V[IJm1K])/(p->DYP[JP]+p->DYP[JM1]))
    
         + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(V[IJKp1] - V[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

double nhflow_gradient::dvdz(double *V)
{
	grad = p->sigz[IJ]*(V[IJKp1] - V[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
    
    if(k==0)
    grad = p->sigz[IJ]*(V[IJK] - V[IJKm1])/(p->DZP[KM1]);
    
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


