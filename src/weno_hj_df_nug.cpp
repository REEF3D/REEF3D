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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"weno_hj_df_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"

weno_hj_df_nug::weno_hj_df_nug(lexer* p):weno_nug_func(p),tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
    if(p->B269==0 && p->D11!=4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269==0 && p->D11==4)
    pflux = new flux_HJ_CDS4(p);
    
    if(p->B269>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans(p);
}

weno_hj_df_nug::~weno_hj_df_nug()
{
}

void weno_hj_df_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;
    
    if(ipol==1)
    {
    uf=1;
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXN,p->DYP,p->DZP);
    }

    if(ipol==2)
    {
    vf=1;
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXP,p->DYN,p->DZP);
    }

    if(ipol==3)
    {
    wf=1;
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXP,p->DYP,p->DZN);
    }

    if(ipol==4)
    FLUIDLOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXP,p->DYP,p->DZP);
    
    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DXP,p->DYP,p->DZP);
}

double weno_hj_df_nug::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DXD,double *DYD, double *DZD)
{
        DX=DXD;
        DY=DYD;
        DZ=DZD;
        
		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
		L = -iadvec*fx(p,a,b,uvel,ipol,iadvec);
        
        if(p->j_dir==1)
        L -= jadvec*fy(p,a,b,vvel,ipol,jadvec);
        
        L -= kadvec*fz(p,a,b,wvel,ipol,kadvec);
        
		return L;
}

double weno_hj_df_nug::fx(lexer *p,fdm *a, field& b, field& uvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
    if(p->X45==0)
	iqmin_0(p,a,b,uvel,ipol);
    
    if(p->X45==1)
	iqmin_1(p,a,b,uvel,ipol);
    
    if(p->X45==2)
	iqmin_2(p,a,b,uvel,ipol);
    
    if(p->X45==3)
	iqmin_3(p,a,b,uvel,ipol);
    
    if(p->X45==4)
	iqmin_4(p,a,b,uvel,ipol);
    
    if(p->X45==5)
	iqmin_5(p,a,b,uvel,ipol);
    
	is_min_x();
	weight_min_x();
    
    
    is();
    alpha();
	weight();
    
    double a1,a2,a3;
            
 
            a1 = cfx[IP][uf][0]/pow(is1x+psi,2.0);
            a2 = cfx[IP][uf][1]/pow(is2x+psi,2.0); 
            a3 = cfx[IP][uf][2]/pow(is3x+psi,2.0); 
    
	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
    if(p->X45==0)
	iqmax_0(p,a,b,uvel,ipol);
    
    if(p->X45==1)
	iqmax_1(p,a,b,uvel,ipol);
    
    if(p->X45==2)
	iqmax_2(p,a,b,uvel,ipol);
    
    if(p->X45==3)
	iqmax_3(p,a,b,uvel,ipol);
    
    if(p->X45==4)
	iqmax_4(p,a,b,uvel,ipol);
    
    if(p->X45==5)
	iqmax_5(p,a,b,uvel,ipol);
    
	is_max_x();
	weight_max_x();
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
	return grad;
}

double weno_hj_df_nug::fy(lexer *p,fdm *a, field& b, field& vvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
    if(p->X45==0)
	jqmin_0(p,a,b,vvel,ipol);
    
    if(p->X45==1)
	jqmin_1(p,a,b,vvel,ipol);
    
    if(p->X45==2)
	jqmin_2(p,a,b,vvel,ipol);
    
    if(p->X45==3)
	jqmin_3(p,a,b,vvel,ipol);
    
    if(p->X45==4)
	jqmin_4(p,a,b,vvel,ipol);
    
    if(p->X45==5)
	jqmin_5(p,a,b,vvel,ipol);
    
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
    if(p->X45==0)
	jqmax_0(p,a,b,vvel,ipol);
    
    if(p->X45==1)
	jqmax_1(p,a,b,vvel,ipol);
    
    if(p->X45==2)
	jqmax_2(p,a,b,vvel,ipol);
    
    if(p->X45==3)
	jqmax_3(p,a,b,vvel,ipol);
    
    if(p->X45==4)
	jqmax_4(p,a,b,vvel,ipol);
    
    if(p->X45==5)
	jqmax_5(p,a,b,vvel,ipol);
    
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}
	
	return grad;
}

double weno_hj_df_nug::fz(lexer *p,fdm *a, field& b, field& wvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
    if(p->X45==0)
	kqmin_0(p,a,b,wvel,ipol);
    
    if(p->X45==1)
	kqmin_1(p,a,b,wvel,ipol);
    
    if(p->X45==2)
	kqmin_2(p,a,b,wvel,ipol);
    
    if(p->X45==3)
	kqmin_3(p,a,b,wvel,ipol);
    
    if(p->X45==4)
	kqmin_4(p,a,b,wvel,ipol);
    
    if(p->X45==5)
	kqmin_5(p,a,b,wvel,ipol);
    
	is_min_z();
	weight_min_z();
	
	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
    if(p->X45==0)
	kqmax_0(p,a,b,wvel,ipol);
    
    if(p->X45==1)
	kqmax_1(p,a,b,wvel,ipol);
    
    if(p->X45==2)
	kqmax_2(p,a,b,wvel,ipol);
    
    if(p->X45==3)
	kqmax_3(p,a,b,wvel,ipol);
    
    if(p->X45==4)
	kqmax_4(p,a,b,wvel,ipol);
    
    if(p->X45==5)
	kqmax_5(p,a,b,wvel,ipol);
    
	is_max_z();
	weight_max_z();
	
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}

	return grad;
}

void weno_hj_df_nug::iqmin_0(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
	q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno_hj_df_nug::jqmin_0(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
	q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno_hj_df_nug::kqmin_0(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
	q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void weno_hj_df_nug::iqmax_0(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
	q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void weno_hj_df_nug::jqmax_0(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
	q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void weno_hj_df_nug::kqmax_0(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
	q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}

void weno_hj_df_nug::iqmin_1(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(fabs(a->fb(i-2,j,k))  > 0.6*p->DXM && fabs(a->fb(i-3,j,k))  > 0.6*p->DXM)
        q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];

    if(fabs(a->fb(i-1,j,k))  > 0.6*p->DXM && fabs(a->fb(i-2,j,k))  > 0.6*p->DXM)
    q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];

    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i-1,j,k))  > 0.6*p->DXM)
        q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];

    if(fabs(a->fb(i+1,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
        q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];

    if(fabs(a->fb(i+2,j,k))  > 0.6*p->DXM && fabs(a->fb(i+1,j,k))  > 0.6*p->DXM)
        q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno_hj_df_nug::jqmin_1(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(fabs(a->fb(i,j-2,k))  > 0.6*p->DXM && fabs(a->fb(i,j-3,k))  > 0.6*p->DXM)
    q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];

    if(fabs(a->fb(i,j-1,k))  > 0.6*p->DXM && fabs(a->fb(i,j-2,k))  > 0.6*p->DXM)
    q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];

    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j-1,k))  > 0.6*p->DXM)
    q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];

    if(fabs(a->fb(i,j+1,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
    q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];

    if(fabs(a->fb(i,j+2,k))  > 0.6*p->DXM && fabs(a->fb(i,j+1,k))  > 0.6*p->DXM)
    q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno_hj_df_nug::kqmin_1(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(fabs(a->fb(i,j,k-2))  > 0.6*p->DXM && fabs(a->fb(i,j,k-3))  > 0.6*p->DXM)
    q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];

    if(fabs(a->fb(i,j,k-1))  > 0.6*p->DXM && fabs(a->fb(i,j,k-2))  > 0.6*p->DXM)
    q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];

    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k-1))  > 0.6*p->DXM)
    q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];

    if(fabs(a->fb(i,j,k+1))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
    q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];

    if(fabs(a->fb(i,j,k+2))  > 0.6*p->DXM && fabs(a->fb(i,j,k+1))  > 0.6*p->DXM)
    q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void weno_hj_df_nug::iqmax_1(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(fabs(a->fb(i-1,j,k))  > 0.6*p->DXM && fabs(a->fb(i-2,j,k))  > 0.6*p->DXM)
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];

    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i-1,j,k))  > 0.6*p->DXM)
    q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];

    if(fabs(a->fb(i+1,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
    q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];

    if(fabs(a->fb(i+2,j,k))  > 0.6*p->DXM && fabs(a->fb(i+1,j,k))  > 0.6*p->DXM)
    q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];

    if(fabs(a->fb(i+3,j,k))  > 0.6*p->DXM && fabs(a->fb(i+2,j,k))  > 0.6*p->DXM)
    q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void weno_hj_df_nug::jqmax_1(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(fabs(a->fb(i,j-1,k))  > 0.6*p->DXM && fabs(a->fb(i,j-2,k))  > 0.6*p->DXM)
    q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];

    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j-1,k))  > 0.6*p->DXM)
    q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];

    if(fabs(a->fb(i,j+1,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
    q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];

    if(fabs(a->fb(i,j+2,k))  > 0.6*p->DXM && fabs(a->fb(i,j+1,k))  > 0.6*p->DXM)
    q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];

    if(fabs(a->fb(i,j+3,k))  > 0.6*p->DXM && fabs(a->fb(i,j+2,k))  > 0.6*p->DXM)
    q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void weno_hj_df_nug::kqmax_1(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(fabs(a->fb(i,j,k-1))  > 0.6*p->DXM && fabs(a->fb(i,j,k-2))  > 0.6*p->DXM)
    q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];

    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k-1))  > 0.6*p->DXM)
    q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];

    if(fabs(a->fb(i,j,k+1))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
    q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];

    if(fabs(a->fb(i,j,k+2))  > 0.6*p->DXM && fabs(a->fb(i,j,k+1))  > 0.6*p->DXM)
    q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];

    if(fabs(a->fb(i,j,k+3))  > 0.6*p->DXM && fabs(a->fb(i,j,k+2))  > 0.6*p->DXM)
    q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}

void weno_hj_df_nug::iqmin_2(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = 0.5*(a->fbh5(i-2,j,k)+a->fbh5(i-3,j,k))*(f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
	q2 = 0.5*(a->fbh5(i-1,j,k)+a->fbh5(i-2,j,k))*(f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q3 = 0.5*(a->fbh5(i,j,k)+a->fbh5(i-1,j,k))*(f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q4 = 0.5*(a->fbh5(i+1,j,k)+a->fbh5(i,j,k))*(f(i+1,j,k)-f(i,j,k))/DX[IP];
	q5 = 0.5*(a->fbh5(i+2,j,k)+a->fbh5(i+1,j,k))*(f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno_hj_df_nug::jqmin_2(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = 0.5*(a->fbh5(i,j-2,k)+a->fbh5(i,j-3,k))*(f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
	q2 = 0.5*(a->fbh5(i,j-1,k)+a->fbh5(i,j-2,k))*(f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q3 = 0.5*(a->fbh5(i,j,k)+a->fbh5(i,j-1,k))*(f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q4 = 0.5*(a->fbh5(i,j+1,k)+a->fbh5(i,j,k))*(f(i,j+1,k)-f(i,j,k))/DY[JP];
	q5 = 0.5*(a->fbh5(i,j+2,k)+a->fbh5(i,j+1,k))*(f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno_hj_df_nug::kqmin_2(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = 0.5*(a->fbh5(i,j,k-2)+a->fbh5(i,j,k-3))*(f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
	q2 = 0.5*(a->fbh5(i,j,k-1)+a->fbh5(i,j,k-2))*(f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q3 = 0.5*(a->fbh5(i,j,k)+a->fbh5(i,j,k-1))*(f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q4 = 0.5*(a->fbh5(i,j,k+1)+a->fbh5(i,j,k))*(f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q5 = 0.5*(a->fbh5(i,j,k+2)+a->fbh5(i,j,k+1))*(f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void weno_hj_df_nug::iqmax_2(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1 = 0.5*(a->fbh5(i-1,j,k)+a->fbh5(i-2,j,k))*(f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = 0.5*(a->fbh5(i,j,k)+a->fbh5(i-1,j,k))*(f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = 0.5*(a->fbh5(i+1,j,k)+a->fbh5(i,j,k))*(f(i+1,j,k)-f(i,j,k))/DX[IP];
	q4 = 0.5*(a->fbh5(i+2,j,k)+a->fbh5(i+1,j,k))*(f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
	q5 = 0.5*(a->fbh5(i+3,j,k)+a->fbh5(i+2,j,k))*(f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void weno_hj_df_nug::jqmax_2(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = 0.5*(a->fbh5(i,j-1,k)+a->fbh5(i,j-2,k))*(f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = 0.5*(a->fbh5(i,j,k)+a->fbh5(i,j-1,k))*(f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = 0.5*(a->fbh5(i,j+1,k)+a->fbh5(i,j,k))*(f(i,j+1,k)-f(i,j,k))/DY[JP];
	q4 = 0.5*(a->fbh5(i,j+2,k)+a->fbh5(i,j+1,k))*(f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
	q5 = 0.5*(a->fbh5(i,j+3,k)+a->fbh5(i,j+2,k))*(f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void weno_hj_df_nug::kqmax_2(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = 0.5*(a->fbh5(i,j,k-1)+a->fbh5(i,j,k-2))*(f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = 0.5*(a->fbh5(i,j,k)+a->fbh5(i,j,k-1))*(f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = 0.5*(a->fbh5(i,j,k+1)+a->fbh5(i,j,k))*(f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q4 = 0.5*(a->fbh5(i,j,k+2)+a->fbh5(i,j,k+1))*(f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
	q5 = 0.5*(a->fbh5(i,j,k+3)+a->fbh5(i,j,k+2))*(f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}


void weno_hj_df_nug::iqmin_3(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(a->fbh4(i-2,j,k)  < 0.5 && a->fbh4(i-3,j,k)  < 0.5)
        q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];

    if(a->fbh4(i-1,j,k)  < 0.5 && a->fbh4(i-2,j,k)  < 0.5)
    q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];

    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i-1,j,k)  < 0.5)
        q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];

    if(a->fbh4(i+1,j,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
        q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];

    if(a->fbh4(i+2,j,k)  < 0.5 && a->fbh4(i+1,j,k)  < 0.5)
        q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno_hj_df_nug::jqmin_3(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(a->fbh4(i,j-2,k)  < 0.5 && a->fbh4(i,j-3,k)  < 0.5)
    q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];

    if(a->fbh4(i,j-1,k)  < 0.5 && a->fbh4(i,j-2,k)  < 0.5)
    q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];

    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j-1,k)  < 0.5)
    q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];

    if(a->fbh4(i,j+1,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];

    if(a->fbh4(i,j+2,k)  < 0.5 && a->fbh4(i,j+1,k)  < 0.5)
    q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno_hj_df_nug::kqmin_3(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(a->fbh4(i,j,k-2)  < 0.5 && a->fbh4(i,j,k-3)  < 0.5)
    q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];

    if(a->fbh4(i,j,k-1)  < 0.5 && a->fbh4(i,j,k-2)  < 0.5)
    q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];

    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j,k-1)  < 0.5)
    q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];

    if(a->fbh4(i,j,k+1)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];

    if(a->fbh4(i,j,k+2)  < 0.5 && a->fbh4(i,j,k+1)  < 0.5)
    q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void weno_hj_df_nug::iqmax_3(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(a->fbh4(i-1,j,k)  < 0.5 && a->fbh4(i-2,j,k)  < 0.5)
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];

    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i-1,j,k)  < 0.5)
    q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];

    if(a->fbh4(i+1,j,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];

    if(a->fbh4(i+2,j,k)  < 0.5 && a->fbh4(i+1,j,k)  < 0.5)
    q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];

    if(a->fbh4(i+3,j,k)  < 0.5 && a->fbh4(i+2,j,k)  < 0.5)
    q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void weno_hj_df_nug::jqmax_3(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(a->fbh4(i,j-1,k)  < 0.5 && a->fbh4(i,j-2,k)  < 0.5)
    q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];

    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j-1,k)  < 0.5)
    q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];

    if(a->fbh4(i,j+1,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];

    if(a->fbh4(i,j+2,k)  < 0.5 && a->fbh4(i,j+1,k)  < 0.5)
    q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];

    if(a->fbh4(i,j+3,k)  < 0.5 && a->fbh4(i,j+2,k)  < 0.5)
    q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void weno_hj_df_nug::kqmax_3(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;

    if(a->fbh4(i,j,k-1)  < 0.5 && a->fbh4(i,j,k-2)  < 0.5)
    q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];

    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j,k-1)  < 0.5)
    q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];

    if(a->fbh4(i,j,k+1)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];

    if(a->fbh4(i,j,k+2)  < 0.5 && a->fbh4(i,j,k+1)  < 0.5)
    q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];

    if(a->fbh4(i,j,k+3)  < 0.5 && a->fbh4(i,j,k+2)  < 0.5)
    q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}

void weno_hj_df_nug::iqmin_4(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = (a->fbh5(i-3,j,k))*(f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
	q2 = (a->fbh5(i-2,j,k))*(f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q3 = (a->fbh5(i-1,j,k))*(f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q4 = (a->fbh5(i,j,k))*(f(i+1,j,k)-f(i,j,k))/DX[IP];
	q5 = (a->fbh5(i+1,j,k))*(f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno_hj_df_nug::jqmin_4(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = (a->fbh5(i,j-3,k))*(f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
	q2 = (a->fbh5(i,j-2,k))*(f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q3 = (a->fbh5(i,j-1,k))*(f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q4 = (a->fbh5(i,j,k))*(f(i,j+1,k)-f(i,j,k))/DY[JP];
	q5 = (a->fbh5(i,j+1,k))*(f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno_hj_df_nug::kqmin_4(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = (a->fbh5(i,j,k-3))*(f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
	q2 = (a->fbh5(i,j,k-2))*(f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q3 = (a->fbh5(i,j,k-1))*(f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q4 = (a->fbh5(i,j,k))*(f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q5 = (a->fbh5(i,j,k+1))*(f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void weno_hj_df_nug::iqmax_4(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1 = (a->fbh5(i-1,j,k))*(f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (a->fbh5(i,j,k))*(f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (a->fbh5(i+1,j,k))*(f(i+1,j,k)-f(i,j,k))/DX[IP];
	q4 = (a->fbh5(i+2,j,k))*(f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
	q5 = (a->fbh5(i+3,j,k))*(f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void weno_hj_df_nug::jqmax_4(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = (a->fbh5(i,j-1,k))*(f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (a->fbh5(i,j,k))*(f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (a->fbh5(i,j+1,k))*(f(i,j+1,k)-f(i,j,k))/DY[JP];
	q4 = (a->fbh5(i,j+2,k))*(f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
	q5 = (a->fbh5(i,j+3,k))*(f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void weno_hj_df_nug::kqmax_4(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = (a->fbh5(i,j,k-1))*(f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (a->fbh5(i,j,k))*(f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (a->fbh5(i,j,k+1))*(f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q4 = (a->fbh5(i,j,k+2))*(f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
	q5 = (a->fbh5(i,j,k+3))*(f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}

void weno_hj_df_nug::iqmin_5(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i-2,j,k)>0.0 && a->topo(i-2,j,k)>0.0 && a->solid(i-3,j,k)>0.0 && a->topo(i-3,j,k)>0.0)
	q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
    
    if(a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0 && a->solid(i-2,j,k)>0.0 && a->topo(i-2,j,k)>0.0)
	q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0)
	q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
    
    if(a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
    
    if(a->solid(i+2,j,k)>0.0 && a->topo(i+2,j,k)>0.0 && a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0)
	q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void weno_hj_df_nug::jqmin_5(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j-2,k)>0.0 && a->topo(i,j-2,k)>0.0 && a->solid(i,j-3,k)>0.0 && a->topo(i,j-3,k)>0.0)
	q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
    
    if(a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0 && a->solid(i,j-2,k)>0.0 && a->topo(i,j-2,k)>0.0)
	q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0)
	q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
    
    if(a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
    
    if(a->solid(i,j+2,k)>0.0 && a->topo(i,j+2,k)>0.0 && a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0)
	q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void weno_hj_df_nug::kqmin_5(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j,k-2)>0.0 && a->topo(i,j,k-2)>0.0 && a->solid(i,j,k-3)>0.0 && a->topo(i,j,k-3)>0.0)
	q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
    
    if(a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0 && a->solid(i,j,k-2)>0.0 && a->topo(i,j,k-2)>0.0)
	q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0)
	q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
    
    if(a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
    
    if(a->solid(i,j,k+2)>0.0 && a->topo(i,j,k+2)>0.0 && a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0)
	q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void weno_hj_df_nug::iqmax_5(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0 && a->solid(i-2,j,k)>0.0 && a->topo(i-2,j,k)>0.0)
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0)
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
    
    if(a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
    
    if(a->solid(i+2,j,k)>0.0 && a->topo(i+2,j,k)>0.0 && a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0)
	q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
    
    if(a->solid(i+3,j,k)>0.0 && a->topo(i+3,j,k)>0.0 && a->solid(i+2,j,k)>0.0 && a->topo(i+2,j,k)>0.0)
	q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void weno_hj_df_nug::jqmax_5(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0 && a->solid(i,j-2,k)>0.0 && a->topo(i,j-2,k)>0.0)
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0)
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
    
    if(a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
    
    if(a->solid(i,j+2,k)>0.0 && a->topo(i,j+2,k)>0.0 && a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0)
	q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
    
    if(a->solid(i,j+3,k)>0.0 && a->topo(i,j+3,k)>0.0 && a->solid(i,j+2,k)>0.0 && a->topo(i,j+2,k)>0.0)
	q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void weno_hj_df_nug::kqmax_5(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0 && a->solid(i,j,k-2)>0.0 && a->topo(i,j,k-2)>0.0)
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0)
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
    
    if(a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
    
    if(a->solid(i,j,k+2)>0.0 && a->topo(i,j,k+2)>0.0 && a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0)
	q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
    
    if(a->solid(i,j,k+3)>0.0 && a->topo(i,j,k+3)>0.0 && a->solid(i,j,k+2)>0.0 && a->topo(i,j,k+2)>0.0)
	q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}

void weno_hj_df_nug::is()
{
	is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
}

void weno_hj_df_nug::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void weno_hj_df_nug::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
