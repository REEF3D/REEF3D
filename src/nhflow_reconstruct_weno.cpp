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

#include"nhflow_reconstruct_weno.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_reconstruct_weno::nhflow_reconstruct_weno(lexer* p, patchBC_interface *ppBC) : weno_nug_func(p), dfdx(p), dfdy(p)
{
    pBC = ppBC;
    p->Darray(DFDX,p->imax*p->jmax*(p->kmax+2));
    
    uf=vf=wf=0;
}

nhflow_reconstruct_weno::~nhflow_reconstruct_weno()
{
}

void nhflow_reconstruct_weno::reconstruct_2D_x(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fs, slice &fn)
{
    SLICELOOP4
    dfdx(i,j) = 0.0;
    
    // gradient
    SLICELOOP4
    WETDRY
    {
    dfdx_plus = (f(i+1,j) - f(i,j))/p->DXP[IP];
    dfdx_min  = (f(i,j) - f(i-1,j))/p->DXP[IM1];
    
    dfdx(i,j) = limiter(dfdx_plus,dfdx_min);
    }
    
    pgc->gcsl_start1(p,dfdx,1);
    
    // reconstruct
    SLICELOOP1  
    {
    fs(i,j) = f(i,j)   + 0.5*p->DXP[IP]*dfdx(i,j); 
    fn(i,j) = f(i+1,j) - 0.5*p->DXP[IP1]*dfdx(i+1,j);
    }
    
    
    SLICELOOP1
    if(p->deep[IJ]==1 && p->deep[Ip1J]==1)
    {
    // left
	iqmin_sl(p,f);
	is_min_x();
	weight_min_x();

	fs(i,j) = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
            + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	
    // right
	iqmax_sl(p,f);
	is_max_x();
	weight_max_x();
    
	fn(i,j) = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
            + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
    pgc->gcsl_start1(p,fs,10);
    pgc->gcsl_start1(p,fn,10);
}

void nhflow_reconstruct_weno::reconstruct_2D_y(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fe, slice &fw)
{
    if(p->j_dir==1)
    {
    SLICELOOP4
    dfdy(i,j) = 0.0;
    
    // gradient
    SLICELOOP4
    WETDRY
    {
    dfdy_plus = (f(i,j+1) - f(i,j))/p->DYP[JP];
    dfdy_min  = (f(i,j) - f(i,j-1))/p->DYP[JM1];
    
    dfdy(i,j) = limiter(dfdy_plus,dfdy_min);
    }
    
    pgc->gcsl_start2(p,dfdy,1);
    
    // reconstruct
    
    SLICELOOP2 
    {
    fe(i,j) = f(i,j)   + 0.5*p->DYP[JP]*dfdy(i,j); 
    fw(i,j) = f(i,j+1) - 0.5*p->DYP[JP1]*dfdy(i,j+1); 
    }
    
    pgc->gcsl_start2(p,fe,1);
    pgc->gcsl_start2(p,fw,1);
    }
    
    
    
    if(p->j_dir==1)
    SLICELOOP2
    if(p->deep[IJ]==1 && p->deep[IJp1]==1)
	{
	jqmin_sl(p,f);
	is_min_y();
	weight_min_y();
	
	fe(i,j) = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
            + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
            + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));

	jqmax_sl(p,f);
	is_max_y();
	weight_max_y();
	
	fw(i,j) = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
            + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
            + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}
    
    pgc->gcsl_start2(p,fe,11);
    pgc->gcsl_start2(p,fw,11);
}

void nhflow_reconstruct_weno::reconstruct_2D_WL(lexer* p, ghostcell *pgc, fdm_nhf *d)
{
    // water level  
    SLICELOOP1
    d->dfx(i,j) = 0.5*(d->depth(i+1,j)+d->depth(i,j));
    
    SLICELOOP2
    d->dfy(i,j) = 0.5*(d->depth(i,j+1)+d->depth(i,j));
    
    pgc->gcsl_start1(p,d->dfx,1);
    pgc->gcsl_start2(p,d->dfy,1);

    
    SLICELOOP1
    {
    d->Ds(i,j) = MAX(d->ETAs(i,j) + 0.5*(d->depth(i+1,j)+d->depth(i,j)), p->A544);
    d->Dn(i,j) = MAX(d->ETAn(i,j) + 0.5*(d->depth(i+1,j)+d->depth(i,j)), p->A544);
    }
    
    if(p->j_dir==1)
    SLICELOOP2
    {
    d->De(i,j) = MAX(d->ETAe(i,j)  + 0.5*(d->depth(i,j+1)+d->depth(i,j)), p->A544);
    d->Dw(i,j) = MAX(d->ETAw(i,j)  + 0.5*(d->depth(i,j+1)+d->depth(i,j)), p->A544);
    }
    
    pgc->gcsl_start1(p,d->Ds,1);
    pgc->gcsl_start1(p,d->Dn,1);
    pgc->gcsl_start2(p,d->De,1);
    pgc->gcsl_start2(p,d->Dw,1);
}

void nhflow_reconstruct_weno::reconstruct_3D_x(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fx, double *Fs, double *Fn)
{
    LOOP
    DFDX[IJK] = 0.0;
    
    LOOP
    WETDRY
    {
    dfdx_plus = (Fx[Ip1JK] - Fx[IJK])/p->DXP[IP];
    dfdx_min  = (Fx[IJK] - Fx[Im1JK])/p->DXP[IM1];
    
    DFDX[IJK] = limiter(dfdx_plus,dfdx_min);
    }

    pgc->start1V(p,DFDX,1);

    // reconstruct
    ULOOP 
    {
    Fs[IJK] = (Fx[IJK]    + 0.5*p->DXP[IP]*DFDX[IJK]); 
    Fn[IJK] = (Fx[Ip1JK]  - 0.5*p->DXP[IP1]*DFDX[Ip1JK]);
    }
    
    
    
    ULOOP
    if(p->deep[IJ]==1 && p->deep[Ip1J]==1)
    {
    // left
	iqmin(p,Fx);
	is_min_x();
	weight_min_x();

	Fs[IJK] =     w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
                + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
                + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	
    // right
	iqmax(p,Fx);
	is_max_x();
	weight_max_x();
    
	Fn[IJK] =     w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
                + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
                + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
    pgc->start1V(p,Fs,10);
    pgc->start1V(p,Fn,10);

}

void nhflow_reconstruct_weno::reconstruct_3D_y(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fy, double *Fe, double *Fw)
{
    if(p->j_dir==1)
    {
    LOOP
    DFDX[IJK] = 0.0;
    
    LOOP
    WETDRY
    {
    dfdy_plus = (Fy[IJp1K] - Fy[IJK])/p->DYP[JP];
    dfdy_min  = (Fy[IJK] - Fy[IJm1K])/p->DYP[JM1];

    DFDX[IJK] = limiter(dfdy_plus,dfdy_min);
    }
    
    pgc->start2V(p,DFDX,1);
    
    // reconstruct
    VLOOP
    {
    Fe[IJK] = (Fy[IJK]    + 0.5*p->DYP[JP]*DFDX[IJK]); 
    Fw[IJK] = (Fy[IJp1K]  - 0.5*p->DYP[JP1]*DFDX[IJp1K]);
    }
    }
    
    if(p->j_dir==1)
    VLOOP
    if(p->deep[IJ]==1 && p->deep[IJp1]==1)
	{
	jqmin(p,Fy);
	is_min_y();
	weight_min_y();
	
	Fe[IJK] =     w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
                + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
                + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));

	jqmax(p,Fy);
	is_max_y();
	weight_max_y();
	
	Fw[IJK] =    w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
                + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
                + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}
    
    pgc->start2V(p,Fe,11);
    pgc->start2V(p,Fw,11);
}

void nhflow_reconstruct_weno::reconstruct_3D_z(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fz, double *Fb, double *Ft)
{
    // gradient
    LOOP
    DFDX[IJK] = 0.0;
    
    LOOP
    WETDRY
    {
    dfdz_plus = (Fz[IJKp1] - Fz[IJK])/p->DZP[KP];
    dfdz_min  = (Fz[IJK] - Fz[IJKm1])/p->DZP[KM1];
    
    DFDX[IJK] = limiter(dfdz_plus,dfdz_min);
    }
    
    pgc->start3V(p,DFDX,1);
    
    // reconstruct
    WLOOP 
    {
    Fb[IJK] = (Fz[IJK]    + 0.5*p->DZN[KP]*DFDX[IJK]); 
    Ft[IJK] = (Fz[IJKp1]  - 0.5*p->DZN[KP1]*DFDX[IJKp1]);
    }
}
/*
void nhflow_reconstruct_weno::reconstruct_3D_z(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fz, double *Fb, double *Ft)
{
    LOOP
    DFDX[IJK] = 0.0;
    
    LOOP
    WETDRY
    {
    dfdx_plus = (Fz[IJKp1] - Fz[IJK])/p->DZP[KP];
    dfdx_min  = (Fz[IJK] - Fz[IJKp1])/p->DZP[KM1];
    
    DFDX[IJK] = limiter(dfdx_plus,dfdx_min);
    }

    pgc->start3V(p,DFDX,1);

    // reconstruct
    WLOOP 
    {
    Fb[IJK] = (Fz[IJK]    + 0.5*p->DZN[KP]*DFDX[IJK]); 
    Ft[IJK] = (Fz[IJKp1]  - 0.5*p->DZN[KP1]*DFDX[IJKp1]);
    }
    
    
    
    WLOOP
    if(p->deep[IJ]==1)
    {
    // left
	kqmin(p,Fz);
	is_min_z();
	weight_min_z();

	Fb[IJK] =     w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
                + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
                + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	
    // right
	kqmax(p,Fz);
	is_max_z();
	weight_max_z();
    
	Ft[IJK] =     w1x*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
                + w2x*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
                + w3x*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}
    
    pgc->start3V(p,Fb,12);
    pgc->start3V(p,Ft,12);
}*/

inline void nhflow_reconstruct_weno::iqmin(lexer *p, double *F)
{	 
    q1 = F[Im2JK];
    q2 = F[Im1JK];
    q3 = F[IJK];
    q4 = F[Ip1JK];
    q5 = F[Ip2JK];
}

inline void nhflow_reconstruct_weno::jqmin(lexer *p, double *F)
{
    q1 = F[IJm2K];
    q2 = F[IJm1K];
    q3 = F[IJK];
    q4 = F[IJp1K];
    q5 = F[IJp2K];
}

inline void nhflow_reconstruct_weno::kqmin(lexer *p, double *F)
{
    q1 = F[IJKm2];
    q2 = F[IJKm1];
    q3 = F[IJK];
    q4 = F[IJKp1];
    q5 = F[IJKp2];
}

inline void nhflow_reconstruct_weno::iqmax(lexer *p, double *F)
{
    q1 = F[Im1JK];
    q2 = F[IJK];
    q3 = F[Ip1JK];
    q4 = F[Ip2JK];
    q5 = F[Ip3JK];
}

inline void nhflow_reconstruct_weno::jqmax(lexer *p, double *F)
{
    q1 = F[IJm1K];
    q2 = F[IJK];
    q3 = F[IJp1K];
    q4 = F[IJp2K];
    q5 = F[IJp3K];
}

inline void nhflow_reconstruct_weno::kqmax(lexer *p, double *F)
{
	q1 = F[IJKm1];
    q2 = F[IJK];
    q3 = F[IJKp1];
    q4 = F[IJKp2];
    q5 = F[IJKp3];
}

inline void nhflow_reconstruct_weno::iqmin_sl(lexer *p, slice& f)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

inline void nhflow_reconstruct_weno::jqmin_sl(lexer *p, slice& f)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

inline void nhflow_reconstruct_weno::iqmax_sl(lexer *p, slice& f)
{
	q1 = f(i-1,j);
	q2 = f(i,j);
	q3 = f(i+1,j);
	q4 = f(i+2,j);
	q5 = f(i+3,j);
}

inline void nhflow_reconstruct_weno::jqmax_sl(lexer *p, slice& f)
{
	q1 = f(i,j-1);
	q2 = f(i,j);
	q3 = f(i,j+1);
	q4 = f(i,j+2);
	q5 = f(i,j+3);
}

inline double nhflow_reconstruct_weno::limiter(double v1, double v2)
{
    val=0.0;
    
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);

    if(r<0.0)
    phi = 0.0;
    
    if(r>=0.0 && r<0.5)
    phi = 2.0*r;
    
    if(r>=0.5 && r<1.0)
    phi = 1.0;
    
    if(r>=1.0)
    phi = MIN(MIN(r,2.0), 2.0/(1.0+r));
    
    val = 0.5*phi*(v1+v2);

    
    return val;
}

