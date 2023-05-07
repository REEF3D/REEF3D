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

#include"nhflow_weno_flux.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_flux_face_cds2.h"

nhflow_weno_flux::nhflow_weno_flux(lexer* p):weno_nug_func(p)
{
        pflux = new nhflow_flux_face_cds2(p);
}

nhflow_weno_flux::~nhflow_weno_flux()
{
}

void nhflow_weno_flux::start(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W)
{
    uf=vf=wf=0;
    
        if(ipol==1)
        LOOP
        d->F[IJK]+=aij(p,d,F,1,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==2)
        LOOP
        d->G[IJK]+=aij(p,d,F,2,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==3)
        LOOP
        d->H[IJK]+=aij(p,d,F,3,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==4)
        LOOP
        d->L[IJK]+=aij(p,d,F,4,U,V,W,p->DXN,p->DYN,p->DZN);
}

double nhflow_weno_flux::aij(lexer* p,fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, double *DX,double *DY, double *DZ)
{
        // convective flux
        pflux->u_flux(d,ipol,U,ivel1,ivel2);
        pflux->v_flux(d,ipol,V,jvel1,jvel2);
        pflux->w_flux(d,ipol,d->omega,kvel1,kvel2);
        

        
        // build stencils
        fv1=fv2=0.0;
		
		i-=1;
		fu1 = fx(p,d,F,U,ipol,ivel1);
		i+=1;
		
		fu2 = fx(p,d,F,U,ipol,ivel2);


		if(p->j_dir==1)
         {
		j-=1;
		fv1 = fy(p,d,F,V,ipol,jvel1);
		j+=1;
		
		fv2 = fy(p,d,F,V,ipol,jvel2);
         }


		k-=1;
		fw1 = fz(p,d,F,W,ipol,kvel1);
		k+=1;
		
		fw2 = fz(p,d,F,W,ipol,kvel2);
    

		L =   - ((ivel2*fu2-ivel1*fu1)/DX[IP]) 
		      - ((jvel2*fv2-jvel1*fv1)/DY[JP]) 
			  - ((kvel2*fw2-kvel1*fw1)/DZ[KP]);
        
		return L;
}

double nhflow_weno_flux::fx(lexer *p, fdm_nhf *d, double *B, double *U, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,B,U,ipol);
	is_min_x();
	weight_min_x();

	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	iqmax(p,B,U,ipol);
	is_max_x();
	weight_max_x();
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
	return grad;
}

double nhflow_weno_flux::fy(lexer *p, fdm_nhf *d, double *B, double *V, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,B,V,ipol);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	jqmax(p,B,V,ipol);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}

	return grad;
}

double nhflow_weno_flux::fz(lexer *p, fdm_nhf *d, double *B, double *W, int ipol, double advec)
{
    grad = 0.0;
    
    double gz1,gz2,gz3;
    double g1,g2,g3;

	if(advec>0.0)
	{
	kqmin(p,B,W,ipol);
	is_min_z();
	weight_min_z();
	
    
	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	kqmax(p,B,W,ipol);
	is_max_z();
	weight_max_z();
	
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}

	return grad;
}

void nhflow_weno_flux::iqmin(lexer *p, double *F, double *U, int ipol)
{	 
    q1 = F[Im2JK];
    q2 = F[Im1JK];
    q3 = F[IJK];
    q4 = F[Ip1JK];
    q5 = F[Ip2JK];
}

void nhflow_weno_flux::jqmin(lexer *p, double *F, double *V, int ipol)
{
    q1 = F[IJm2K];
    q2 = F[IJm1K];
    q3 = F[IJK];
    q4 = F[IJp1K];
    q5 = F[IJp2K];
}

void nhflow_weno_flux::kqmin(lexer *p, double *F, double *W, int ipol)
{
    q1 = F[IJKm2];
    q2 = F[IJKm1];
    q3 = F[IJK];
    q4 = F[IJKp1];
    q5 = F[IJKp2];
}

void nhflow_weno_flux::iqmax(lexer *p, double *F, double *U, int ipol)
{
    q1 = F[Im1JK];
    q2 = F[IJK];
    q3 = F[Ip1JK];
    q4 = F[Ip2JK];
    q5 = F[Ip3JK];
}

void nhflow_weno_flux::jqmax(lexer *p, double *F, double *V, int ipol)
{
    q1 = F[IJm1K];
    q2 = F[IJK];
    q3 = F[IJp1K];
    q4 = F[IJp2K];
    q5 = F[IJp3K];
}

void nhflow_weno_flux::kqmax(lexer *p, double *F, double *W, int ipol)
{
	q1 = F[IJKm1];
    q2 = F[IJK];
    q3 = F[IJKp1];
    q4 = F[IJKp2];
    q5 = F[IJKp3];
}


