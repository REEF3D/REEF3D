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

#include"nhflow_fsf_reconstruct_WENO.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_fsf_reconstruct_weno::nhflow_fsf_reconstruct_weno(lexer* p, patchBC_interface *ppBC) : weno_nug_func(p), dfdx(p), dfdy(p)
{
    pBC = ppBC;
    
    p->Darray(DFDX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DFDY,p->imax*p->jmax*(p->kmax+2));

}

nhflow_fsf_reconstruct_weno::~nhflow_fsf_reconstruct_weno()
{
}

void nhflow_fsf_reconstruct_weno::reconstruct_2D(lexer* p, ghostcell *pgc, fdm_nhf *d, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
{
    ULOOP
    {
    // left
	iqmin_sl(p,f);
	is_min_x();
	weight_min_x();

	fs(i,j) = (w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
            + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2)));
	
    // right
	iqmax_sl(p,f);
	is_max_x();
	weight_max_x();
    
	fn(i,j) = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
            + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
    
    VLOOP
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
    
    pgc->gcsl_start1(p,fs,10);
    pgc->gcsl_start1(p,fn,10);
    pgc->gcsl_start2(p,fe,11);
    pgc->gcsl_start2(p,fw,11);
}

void nhflow_fsf_reconstruct_weno::reconstruct_3D(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fx, double *Fy, double *Fs, double *Fn, double *Fe, double *Fw)
{
    ULOOP
    {
    // left
	iqmin(p,Fx);
	is_min_x();
	weight_min_x();

	Fs[IJK] =        (w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
                    + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
                    + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2)));
	
    // right
	iqmax(p,Fx);
	is_max_x();
	weight_max_x();
    
	Fn[IJK] =            (w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
                        + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
                        + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2)));
	}
    
    
    VLOOP
	{
	jqmin(p,Fy);
	is_min_y();
	weight_min_y();
	
	Fe[IJK] =           (w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
                        + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
                        + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2)));

	jqmax(p,Fy);
	is_max_y();
	weight_max_y();
	
	Fw[IJK] =           (w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
                        + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
                        + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2)));
	}
    
    
    
	pgc->start1V(p,Fs,10);
    pgc->start1V(p,Fs,10);
    pgc->start2V(p,Fe,11);
    pgc->start2V(p,Fw,11);
}

double nhflow_fsf_reconstruct_weno::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}

void nhflow_fsf_reconstruct_weno::iqmin(lexer *p, double *F)
{	 
    q1 = F[Im2JK];
    q2 = F[Im1JK];
    q3 = F[IJK];
    q4 = F[Ip1JK];
    q5 = F[Ip2JK];
}

void nhflow_fsf_reconstruct_weno::jqmin(lexer *p, double *F)
{
    q1 = F[IJm2K];
    q2 = F[IJm1K];
    q3 = F[IJK];
    q4 = F[IJp1K];
    q5 = F[IJp2K];
}

void nhflow_fsf_reconstruct_weno::kqmin(lexer *p, double *F)
{
    q1 = F[IJKm2];
    q2 = F[IJKm1];
    q3 = F[IJK];
    q4 = F[IJKp1];
    q5 = F[IJKp2];
}

void nhflow_fsf_reconstruct_weno::iqmax(lexer *p, double *F)
{
    q1 = F[Im1JK];
    q2 = F[IJK];
    q3 = F[Ip1JK];
    q4 = F[Ip2JK];
    q5 = F[Ip3JK];
}

void nhflow_fsf_reconstruct_weno::jqmax(lexer *p, double *F)
{
    q1 = F[IJm1K];
    q2 = F[IJK];
    q3 = F[IJp1K];
    q4 = F[IJp2K];
    q5 = F[IJp3K];
}

void nhflow_fsf_reconstruct_weno::kqmax(lexer *p, double *F)
{
	q1 = F[IJKm1];
    q2 = F[IJK];
    q3 = F[IJKp1];
    q4 = F[IJKp2];
    q5 = F[IJKp3];
}

void nhflow_fsf_reconstruct_weno::iqmin_sl(lexer *p, slice& f)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

void nhflow_fsf_reconstruct_weno::jqmin_sl(lexer *p, slice& f)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

void nhflow_fsf_reconstruct_weno::iqmax_sl(lexer *p, slice& f)
{
	q1 = f(i+3,j);
	q2 = f(i+2,j);
	q3 = f(i+1,j);
	q4 = f(i,j);
	q5 = f(i-1,j);
}

void nhflow_fsf_reconstruct_weno::jqmax_sl(lexer *p, slice& f)
{
	q1 = f(i,j+3);
	q2 = f(i,j+2);
	q3 = f(i,j+1);
	q4 = f(i,j);
	q5 = f(i,j-1);
}

