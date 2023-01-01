/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"weno3_hj.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"

weno3_hj::weno3_hj(lexer* p) : weno3_nug_func(p)
{
    if(p->B269==0 && p->D11!=4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269==0 && p->D11==4)
    pflux = new flux_HJ_CDS4(p);
    
    if(p->B269>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans(p);
}

weno3_hj::~weno3_hj()
{
}

void weno3_hj::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;
    
    if(ipol==1)
    {
    uf=1;
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP,p->DYN,p->DZN);
    }

    if(ipol==2)
    {
    vf=1;
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN,p->DYP,p->DZN);
    }

    if(ipol==3)
    {
    wf=1;
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN,p->DYN,p->DZP);
    }

    if(ipol==4)
    FLUIDLOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);
    
    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);
}

double weno3_hj::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DXD,double *DYD, double *DZD)
{
        DX=DXD;
        DY=DYD;
        DZ=DZD;
        
		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
		L = -iadvec*fx(p,a,b,uvel,ipol,iadvec) - jadvec*fy(p,a,b,vvel,ipol,jadvec) - kadvec*fz(p,a,b,wvel,ipol,kadvec);
        
		return L;
}

double weno3_hj::fx(lexer *p,fdm *a, field& b, field& uvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,a,b,uvel,ipol);
	is_min_x();
	weight_min_x();

	grad = w1x*(qfx[IP][uf][0][0]*q2 + qfx[IP][uf][0][1]*q3)
    
         + w2x*(qfx[IP][uf][1][0]*q2 - qfx[IP][uf][1][1]*q1);
	}

	if(advec<0.0)
	{
	iqmax(p,a,b,uvel,ipol);
	is_max_x();
	weight_max_x();
    
	grad = w1x*(qfx[IP][uf][2][0]*q2 - qfx[IP][uf][2][1]*q3)
    
         + w2x*(qfx[IP][uf][3][0]*q1 + qfx[IP][uf][3][1]*q2);
	}
    
	return grad;
}

double weno3_hj::fy(lexer *p,fdm *a, field& b, field& vvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,a,b,vvel,ipol);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(qfy[JP][vf][0][0]*q2 + qfy[JP][vf][0][1]*q3)
    
         + w2y*(qfy[JP][vf][1][0]*q2 - qfy[JP][vf][1][1]*q1);
	}

	if(advec<0.0)
	{
	jqmax(p,a,b,vvel,ipol);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(qfy[JP][vf][2][0]*q2 - qfy[JP][vf][2][1]*q3)
    
         + w2y*(qfy[JP][vf][3][0]*q1 + qfy[JP][vf][3][1]*q2);
	}
	
	return grad;
}

double weno3_hj::fz(lexer *p,fdm *a, field& b, field& wvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	kqmin(p,a,b,wvel,ipol);
	is_min_z();
	weight_min_z();
	
	grad = w1z*(qfz[KP][wf][0][0]*q2 + qfz[KP][wf][0][1]*q3)
    
         + w2z*(qfz[KP][wf][1][0]*q2 - qfz[KP][wf][1][1]*q1);
	}

	if(advec<0.0)
	{
	kqmax(p,a,b,wvel,ipol);
	is_max_z();
	weight_max_z();
	
	grad = w1z*(qfz[KP][wf][2][0]*q2 - qfz[KP][wf][2][1]*q3)
    
         + w2z*(qfz[KP][wf][3][0]*q1 + qfz[KP][wf][3][1]*q2);
	}

	return grad;
}

void weno3_hj::iqmin(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
}

void weno3_hj::jqmin(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
}

void weno3_hj::kqmin(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
}

void weno3_hj::iqmax(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
	q1 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q2 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q3 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno3_hj::jqmax(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q2 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q3 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno3_hj::kqmax(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q2 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q3 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

