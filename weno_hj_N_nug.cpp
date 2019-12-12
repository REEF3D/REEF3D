/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"weno_hj_N_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"

weno_hj_N_nug::weno_hj_N_nug(lexer* p):weno_nug_func(p)
{
    if(p->B269==0 && p->D11!=4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269==0 && p->D11==4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans(p);
}

weno_hj_N_nug::~weno_hj_N_nug()
{
}

void weno_hj_N_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;
    
    if(ipol==1)
    {
    fillxvec1(p,a,b);
    uf=1;
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP,p->DYN,p->DZN,a->C1);
    }

    if(ipol==2)
    {
    fillxvec2(p,a,b);
    vf=1;
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN,p->DYP,p->DZN,a->C2);
    }

    if(ipol==3)
    {
    fillxvec3(p,a,b);
    wf=1;
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN,p->DYN,p->DZP,a->C3);
    }

    if(ipol==4)
    {
    fillxvec4(p,a,b);
    FLUIDLOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN,a->C4);
    }
    
    if(ipol==5)
    {
    fillxvec4(p,a,b);
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN,a->C4);
    }
}

double weno_hj_N_nug::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DXD,double *DYD, double *DZD, cpt &C)
{
        DX=DXD;
        DY=DYD;
        DZ=DZD;
        
		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
		L = -iadvec*fx(p,a,b,uvel,ipol,iadvec,C);
        
        if(p->j_dir==1)
        L -= jadvec*fy(p,a,b,vvel,ipol,jadvec,C);
        
        L -= kadvec*fz(p,a,b,wvel,ipol,kadvec,C);
        
		return L;
}

double weno_hj_N_nug::fx(lexer *p,fdm *a, field& b, field& uvel, int ipol, double advec, cpt &C)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,a,b,uvel,C);
	is_min_x();
	weight_min_x();

	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	iqmax(p,a,b,uvel,C);
	is_max_x();
	weight_max_x();
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
	return grad;
}

double weno_hj_N_nug::fy(lexer *p,fdm *a, field& b, field& vvel, int ipol, double advec, cpt &C)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,a,b,vvel,C);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	jqmax(p,a,b,vvel,C);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}
	
	return grad;
}

double weno_hj_N_nug::fz(lexer *p,fdm *a, field& b, field& wvel, int ipol, double advec, cpt &C)
{
    grad = 0.0;

	if(advec>0.0)
	{
	kqmin(p,a,b,wvel,C);
	is_min_z();
	weight_min_z();
	
	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	kqmax(p,a,b,wvel,C);
	is_max_z();
	weight_max_z();
	
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}

	return grad;
}

void weno_hj_N_nug::iqmin(lexer *p,fdm *a, field& f, field& uvel, cpt &C)
{	
	q1 = (a->xvec.V[Im2_J_K] - a->xvec.V[Im3_J_K])/DX[IM3];
	q2 = (a->xvec.V[Im1_J_K] - a->xvec.V[Im2_J_K])/DX[IM2];
	q3 = (a->xvec.V[I_J_K]   - a->xvec.V[Im1_J_K])/DX[IM1];
	q4 = (a->xvec.V[Ip1_J_K] - a->xvec.V[I_J_K]  )/DX[IP];
	q5 = (a->xvec.V[Ip2_J_K] - a->xvec.V[Ip1_J_K])/DX[IP1];
}

void weno_hj_N_nug::jqmin(lexer *p,fdm *a, field& f, field& vvel, cpt &C)
{
	q1 = (a->xvec.V[I_Jm2_K] - a->xvec.V[I_Jm3_K])/DY[JM3];
	q2 = (a->xvec.V[I_Jm1_K] - a->xvec.V[I_Jm2_K])/DY[JM2];
	q3 = (a->xvec.V[I_J_K]   - a->xvec.V[I_Jm1_K])/DY[JM1];
	q4 = (a->xvec.V[I_Jp1_K] - a->xvec.V[I_J_K]  )/DY[JP];
	q5 = (a->xvec.V[I_Jp2_K] - a->xvec.V[I_Jp1_K])/DY[JP1];
}

void weno_hj_N_nug::kqmin(lexer *p,fdm *a, field& f, field& wvel, cpt &C)
{
	q1 = (a->xvec.V[I_J_Km2] - a->xvec.V[I_J_Km3])/DZ[KM3];
	q2 = (a->xvec.V[I_J_Km1] - a->xvec.V[I_J_Km2])/DZ[KM2];
	q3 = (a->xvec.V[I_J_K]   - a->xvec.V[I_J_Km1])/DZ[KM1];
	q4 = (a->xvec.V[I_J_Kp1] - a->xvec.V[I_J_K]  )/DZ[KP];
	q5 = (a->xvec.V[I_J_Kp2] - a->xvec.V[I_J_Kp1])/DZ[KP1];
}

void weno_hj_N_nug::iqmax(lexer *p,fdm *a, field& f, field& uvel, cpt &C)
{
    q1 = (a->xvec.V[Ip3_J_K] - a->xvec.V[Ip2_J_K])/DX[IM2];
	q2 = (a->xvec.V[Ip2_J_K] - a->xvec.V[Ip1_J_K])/DX[IM1];
	q3 = (a->xvec.V[Ip1_J_K] - a->xvec.V[I_J_K]  )/DX[IP];
	q4 = (a->xvec.V[I_J_K]   - a->xvec.V[Im1_J_K])/DX[IP1];
	q5 = (a->xvec.V[Im1_J_K] - a->xvec.V[Im2_J_K])/DX[IP2];
}

void weno_hj_N_nug::jqmax(lexer *p,fdm *a, field& f, field& vvel, cpt &C)
{
	q1 = (a->xvec.V[I_Jp3_K] - a->xvec.V[I_Jp2_K])/DY[JM2];
	q2 = (a->xvec.V[I_Jp2_K] - a->xvec.V[I_Jp1_K])/DY[JM1];
	q3 = (a->xvec.V[I_Jp1_K] - a->xvec.V[I_J_K]  )/DY[JP];
	q4 = (a->xvec.V[I_J_K]   - a->xvec.V[I_Jm1_K])/DY[JP1];
	q5 = (a->xvec.V[I_Jm1_K] - a->xvec.V[I_Jm2_K])/DY[JP2];
}

void weno_hj_N_nug::kqmax(lexer *p,fdm *a, field& f, field& wvel, cpt &C)
{
	q1 = (a->xvec.V[I_J_Kp3] - a->xvec.V[I_J_Kp2])/DZ[KM2];
	q2 = (a->xvec.V[I_J_Kp2] - a->xvec.V[I_J_Kp1])/DZ[KM1];
	q3 = (a->xvec.V[I_J_Kp1] - a->xvec.V[I_J_K]  )/DZ[KP];
	q4 = (a->xvec.V[I_J_K]   - a->xvec.V[I_J_Km1])/DZ[KP1];
	q5 = (a->xvec.V[I_J_Km1] - a->xvec.V[I_J_Km2])/DZ[KP2];
}

