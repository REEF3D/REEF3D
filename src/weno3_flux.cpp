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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"weno3_flux.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

weno3_flux::weno3_flux(lexer* p) : weno3_nug_func(p)
{
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
        
        if(p->D11==3)
        pflux = new flux_face_QOU(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS4(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
        
        if(p->D11==3)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS2(p);
    }
}

weno3_flux::~weno3_flux()
{
}

void weno3_flux::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

double weno3_flux::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{
        pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
		
		i-=1;
		fu1 = fx(p,a,b,uvel,ipol,ivel1);
		i+=1;
		
		fu2 = fx(p,a,b,uvel,ipol,ivel2);

        fv1=fv2=0.0;
		if(p->j_dir==1)
        {
		j-=1;
		fv1 = fy(p,a,b,vvel,ipol,jvel1);
		j+=1;
		
		fv2 = fy(p,a,b,vvel,ipol,jvel2);
        }



		k-=1;
		fw1 = fz(p,a,b,wvel,ipol,kvel1);
		k+=1;
		
		fw2 = fz(p,a,b,wvel,ipol,kvel2);
		
		
		L =   - ((ivel2*fu2-ivel1*fu1)/DX[IP]) 
		      - ((jvel2*fv2-jvel1*fv1)/DY[JP]) 
			  - ((kvel2*fw2-kvel1*fw1)/DZ[KP]);
			  
        
		return L;
}

double weno3_flux::fx(lexer *p,fdm *a, field& b, field& uvel, int ipol, double advec)
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

double weno3_flux::fy(lexer *p,fdm *a, field& b, field& vvel, int ipol, double advec)
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

double weno3_flux::fz(lexer *p,fdm *a, field& b, field& wvel, int ipol, double advec)
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

void weno3_flux::iqmin(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = f(i-1,j,k);
	q2 = f(i,j,k);
	q3 = f(i+1,j,k);
}

void weno3_flux::jqmin(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-1,k);
	q2 = f(i,j,k);
	q3 = f(i,j+1,k);
}

void weno3_flux::kqmin(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-1);
	q2 = f(i,j,k);
	q3 = f(i,j,k+1);
}

void weno3_flux::iqmax(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1 = f(i,j,k);
	q2 = f(i+1,j,k);
	q3 = f(i+2,j,k);
}

void weno3_flux::jqmax(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j,k);
	q2 = f(i,j+1,k);
	q3 = f(i,j+2,k);
}

void weno3_flux::kqmax(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k);
	q2 = f(i,j,k+1);
	q3 = f(i,j,k+2);
}





