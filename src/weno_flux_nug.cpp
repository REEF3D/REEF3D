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

#include"weno_flux_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

weno_flux_nug::weno_flux_nug(lexer* p):weno_nug_func(p)
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

weno_flux_nug::~weno_flux_nug()
{
}

void weno_flux_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

double weno_flux_nug::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{
        pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
        
        fv1=fv2=0.0;
		
		i-=1;
		fu1 = fx(p,a,b,uvel,ipol,ivel1);
		i+=1;
		
		fu2 = fx(p,a,b,uvel,ipol,ivel2);


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

double weno_flux_nug::fx(lexer *p,fdm *a, field& b, field& uvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,b,uvel,ipol);
	is_min_x();
	weight_min_x();

	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	iqmax(p,b,uvel,ipol);
	is_max_x();
	weight_max_x();
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
	return grad;
}

double weno_flux_nug::fy(lexer *p,fdm *a, field& b, field& vvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,b,vvel,ipol);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	jqmax(p,b,vvel,ipol);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}

	return grad;
}

double weno_flux_nug::fz(lexer *p,fdm *a, field& b, field& wvel, int ipol, double advec)
{
    grad = 0.0;
    
    double gz1,gz2,gz3;
    double g1,g2,g3;

	if(advec>0.0)
	{
	kqmin(p,b,wvel,ipol);
	is_min_z();
	weight_min_z();
	
    
	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	kqmax(p,b,wvel,ipol);
	is_max_z();
	weight_max_z();
	
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}

    

	return grad;
}

void weno_flux_nug::iqmin(lexer *p, field& f, field& uvel, int ipol)
{	
	q1 = f(i-2,j,k);
	q2 = f(i-1,j,k);
	q3 = f(i,j,k);
	q4 = f(i+1,j,k);
	q5 = f(i+2,j,k);
}

void weno_flux_nug::jqmin(lexer *p, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-2,k);
	q2 = f(i,j-1,k);
	q3 = f(i,j,k);
	q4 = f(i,j+1,k);
	q5 = f(i,j+2,k);
}

void weno_flux_nug::kqmin(lexer *p, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-2);
	q2 = f(i,j,k-1);
	q3 = f(i,j,k);
	q4 = f(i,j,k+1);
	q5 = f(i,j,k+2);
}

void weno_flux_nug::iqmax(lexer *p, field& f, field& uvel, int ipol)
{
    q1 = f(i-1,j,k);
	q2 = f(i,j,k);
	q3 = f(i+1,j,k);
	q4 = f(i+2,j,k);
	q5 = f(i+3,j,k);
}

void weno_flux_nug::jqmax(lexer *p, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-1,k);
	q2 = f(i,j,k);
	q3 = f(i,j+1,k);
	q4 = f(i,j+2,k);
	q5 = f(i,j+3,k);
}

void weno_flux_nug::kqmax(lexer *p, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-1);
	q2 = f(i,j,k);
	q3 = f(i,j,k+1);
	q4 = f(i,j,k+2);
	q5 = f(i,j,k+3);
}


