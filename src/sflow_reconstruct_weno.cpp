/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sflow_reconstruct_weno.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"patchBC_interface.h"

sflow_reconstruct_weno::sflow_reconstruct_weno(lexer* p, patchBC_interface *ppBC) : weno_nug_func(p), dfdx(p), dfdy(p)
{
    pBC = ppBC;

    uf=vf=wf=0;
}

sflow_reconstruct_weno::~sflow_reconstruct_weno()
{
}

void sflow_reconstruct_weno::reconstruct_x(lexer* p, ghostcell *pgc, fdm2D*, slice& f, slice &fs, slice &fn)
{
    uf=1;
    vf=0;
    wf=0;
    
    SLICELOOP1
    {
    // left
	iqmin(p,f);
	is_min_x();
	weight_min_x();

	fs(i,j) = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
            + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	
    // right
	iqmax(p,f);
	is_max_x();
	weight_max_x();
    
	fn(i,j) = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
            + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
    pgc->gcsl_start1(p,fs,1);
    pgc->gcsl_start1(p,fn,1);
}

void sflow_reconstruct_weno::reconstruct_y(lexer* p, ghostcell *pgc, fdm2D*, slice& f, slice &fe, slice &fw)
{
    uf=0;
    vf=1;
    wf=0;
    
    if(p->j_dir==1)
    SLICELOOP2
	{
	jqmin(p,f);
	is_min_y();
	weight_min_y();
	
	fe(i,j) = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
            + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
            + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));

	jqmax(p,f);
	is_max_y();
	weight_max_y();
	
	fw(i,j) = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
            + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
            + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}
    
    pgc->gcsl_start2(p,fe,1);
    pgc->gcsl_start2(p,fw,1);
}

void sflow_reconstruct_weno::reconstruct_WL(lexer* p, ghostcell *pgc, fdm2D *b)
{
    // water level  
    SLICELOOP1
    b->dfx(i,j) = 0.5*(b->depth(i+1,j)+b->depth(i,j));
    
    SLICELOOP2
    b->dfy(i,j) = 0.5*(b->depth(i,j+1)+b->depth(i,j));
    
    pgc->gcsl_start1(p,b->dfx,1);
    pgc->gcsl_start2(p,b->dfy,1);

    
    SLICELOOP1
    {
    b->Ds(i,j) = MAX(b->ETAs(i,j) + 0.5*(b->depth(i+1,j)+b->depth(i,j)), p->A544);
    b->Dn(i,j) = MAX(b->ETAn(i,j) + 0.5*(b->depth(i+1,j)+b->depth(i,j)), p->A544);
    }
    
    SLICELOOP2
    {
    b->De(i,j) = MAX(b->ETAe(i,j)  + 0.5*(b->depth(i,j+1)+b->depth(i,j)), p->A544);
    b->Dw(i,j) = MAX(b->ETAw(i,j)  + 0.5*(b->depth(i,j+1)+b->depth(i,j)), p->A544);
    }
    
    pgc->gcsl_start1(p,b->Ds,1);
    pgc->gcsl_start1(p,b->Dn,1);
    pgc->gcsl_start2(p,b->De,1);
    pgc->gcsl_start2(p,b->Dw,1);
}

inline void sflow_reconstruct_weno::iqmin(lexer *p, slice& f)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

inline void sflow_reconstruct_weno::jqmin(lexer *p, slice& f)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

inline void sflow_reconstruct_weno::iqmax(lexer *p, slice& f)
{
	q1 = f(i-1,j);
	q2 = f(i,j);
	q3 = f(i+1,j);
	q4 = f(i+2,j);
	q5 = f(i+3,j);
}

inline void sflow_reconstruct_weno::jqmax(lexer *p, slice& f)
{
	q1 = f(i,j-1);
	q2 = f(i,j);
	q3 = f(i,j+1);
	q4 = f(i,j+2);
	q5 = f(i,j+3);
}

