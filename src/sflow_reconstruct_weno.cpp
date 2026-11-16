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

#include"sflow_reconstruct_weno.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"patchBC_interface.h"

sflow_reconstruct_weno::sflow_reconstruct_weno(lexer* p, patchBC_interface *ppBC) : weno_nug_func(p)
{
    pBC = ppBC;
    drs = false;

    uf=vf=wf=0;

    // nonlinear weights (A 213: 0 WENO-JS, 1 WENO-Z, 2 TENO5)
    set_weno_weights(p->A213,p->A222);
}

sflow_reconstruct_weno::sflow_reconstruct_weno(lexer* p, patchBC_interface *ppBC, int own) : weno_nug_func(p,own)
{
    pBC = ppBC;
    drs = false;

    uf=vf=wf=0;

    // nonlinear weights (A 213: 0 WENO-JS, 1 WENO-Z, 2 TENO5)
    set_weno_weights(p->A213,p->A222);
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
    // dry cells in the stencils of the two cells of the face (prepare)
    const int cl = drs ? drx[IJ] : 0;
    const int cr = drs ? drx[IJ + p->jmax] : 0;
    
    // left
	iqmin(p,f);
    
    if(cl==1)
    fs(i,j) = q3;
    
    else
    if(cl==2)
    fs(i,j) = q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3);
    
    else
    {
	is_min_x();
	weight_min_x();

	fs(i,j) = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
            + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
    }
	
    // right
	iqmax(p,f);
    
    if(cr==1)
    fn(i,j) = q3;
    
    else
    if(cr==2)
    fn(i,j) = q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3);
    
    else
    {
	is_max_x();
	weight_max_x();
    
	fn(i,j) = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
            + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
            + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
    }
	}
}

void sflow_reconstruct_weno::reconstruct_y(lexer* p, ghostcell *pgc, fdm2D*, slice& f, slice &fe, slice &fw)
{
    uf=0;
    vf=1;
    wf=0;
    
    if(p->j_dir==1)
    SLICELOOP2
	{
    const int cl = drs ? dry[IJ] : 0;
    const int cr = drs ? dry[IJ + 1] : 0;
    
	jqmin(p,f);
    
    if(cl==1)
    fe(i,j) = q3;
    
    else
    if(cl==2)
    fe(i,j) = q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3);
    
    else
    {
	is_min_y();
	weight_min_y();
	
	fe(i,j) = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
            + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
            + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
    }

	jqmax(p,f);
    
    if(cr==1)
    fw(i,j) = q3;
    
    else
    if(cr==2)
    fw(i,j) = q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3);
    
    else
    {
	is_max_y();
	weight_max_y();
	
	fw(i,j) = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
            + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
            + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
    }
	}
}

// Dry cells in the WENO stencil of a wet cell (wetting and drying A 243 > 0), per direction.
// The values of a dry cell are no data of the flow: its eta is its bed level plus the film A 244, its velocity zero.
// Next to a shoreline the nonlinear weights give the stencils across it no weight, as long as one of the three
// candidate stencils of the cell is wet. Two cases are left, a channel or a passage three cells wide or narrower
// between dry cells:
//  1: every candidate stencil holds a dry cell (dry cells within two cells on both sides, one of them the neighbour):
//     the cell value (first order)
//  2: both outer stencils hold a dry cell and the central one is wet (the middle cell of a channel three cells wide):
//     the central stencil alone. WENO-Z (A 213 1, the default) measures the smoothness with |IS1 - IS3|, which vanishes
//     for the same jump on both sides: the linear fifth-order stencil then ran across both banks (a bank 6 m above the
//     water: 0.1 m in the face levels; the channel had lateral level differences of 0.3 m and no steady state)
// All other cells are unchanged. Solids hold the wet flag of their fluid neighbour (gcsl_start4Vint) and mirrored
// values: not dry.
void sflow_reconstruct_weno::prepare(lexer *p)
{
    drs = false;
    
    if(p->A243==0)
    return;
    
    const int ni = p->imax, nj = p->jmax;
    
    if(int(drx.size())!=ni*nj)
    {
    drx.assign(ni*nj,0);
    dry.assign(ni*nj,0);
    }
    
    auto isdry = [&](int q) { return p->wet[q]==0 && p->flagslice4[q]>0; };
    
    int any = 0;
    
    for(int a=2; a<ni-2; ++a)
    for(int c=0; c<nj; ++c)
    {
    const int q = a*nj + c;
    unsigned char code = 0;
    
        if(p->wet[q]==1)
        {
        const bool m1 = isdry(q-nj), m2 = isdry(q-2*nj), p1 = isdry(q+nj), p2 = isdry(q+2*nj);
        
        // candidates (c-2,c-1,c), (c-1,c,c+1), (c,c+1,c+2)
        const bool d3 = m1 || m2, d2 = m1 || p1, d1 = p1 || p2;
        
        if(d1 && d2 && d3)
        code = 1;
        
        else
        if(d1 && d3)
        code = 2;
        }
    
    drx[q] = code;
    any |= code;
    }
    
    if(p->j_dir==1)
    for(int a=0; a<ni; ++a)
    for(int c=2; c<nj-2; ++c)
    {
    const int q = a*nj + c;
    unsigned char code = 0;
    
        if(p->wet[q]==1)
        {
        const bool m1 = isdry(q-1), m2 = isdry(q-2), p1 = isdry(q+1), p2 = isdry(q+2);
        const bool d3 = m1 || m2, d2 = m1 || p1, d1 = p1 || p2;
        
        if(d1 && d2 && d3)
        code = 1;
        
        else
        if(d1 && d3)
        code = 2;
        }
    
    dry[q] = code;
    any |= code;
    }
    
    drs = any!=0;
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
    b->Ds(i,j) = MAX(b->ETAs(i,j) + b->dfx(i,j), p->A244);
    b->Dn(i,j) = MAX(b->ETAn(i,j) + b->dfx(i,j), p->A244);
    }
    
    SLICELOOP2
    {
    b->De(i,j) = MAX(b->ETAe(i,j) + b->dfy(i,j), p->A244);
    b->Dw(i,j) = MAX(b->ETAw(i,j) + b->dfy(i,j), p->A244);
    }
}

void sflow_reconstruct_weno::iqmin(lexer *p, slice& f)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

void sflow_reconstruct_weno::jqmin(lexer *p, slice& f)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

void sflow_reconstruct_weno::iqmax(lexer *p, slice& f)
{
	q1 = f(i-1,j);
	q2 = f(i,j);
	q3 = f(i+1,j);
	q4 = f(i+2,j);
	q5 = f(i+3,j);
}

void sflow_reconstruct_weno::jqmax(lexer *p, slice& f)
{
	q1 = f(i,j-1);
	q2 = f(i,j);
	q3 = f(i,j+1);
	q4 = f(i,j+2);
	q5 = f(i,j+3);
}

