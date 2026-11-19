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

// The WENO5 reconstruction of one face value from the five cell values q1..q5: the smoothness indicators
// (weno_nug_func::is_min_x), the nonlinear weights (weight_min_x; wt: 0 WENO-JS, 1 WENO-Z, 2 TENO5) and the
// three candidate values, with the same operations in the same order. is: 3 x 3, cf: 3, qf: 3 x 2 coefficients
// of the face; sg = -1 for the left state of a face (the stencil of cell i), +1 for the right state (cell i+1):
// the two differ in the signs of two terms.
template<int WT, int SG>
static inline __attribute__((always_inline)) double sflow_weno5(const double q1, const double q2, const double q3,
                               const double q4, const double q5,
                               const double i00, const double i01, const double i02,
                               const double i10, const double i11, const double i12,
                               const double i20, const double i21, const double i22,
                               const double c1, const double c2, const double c3,
                               const double f00, const double f01, const double f10, const double f11,
                               const double f20, const double f21,
                               const double psi, const double epsilon, const double ct)
{
    const double dq12 = q1 - q2;
    const double dq23 = q2 - q3;
    const double dq32 = q3 - q2;
    const double dq34 = q3 - q4;
    const double dq43 = q4 - q3;
    const double dq54 = q5 - q4;

    const double is1 = i00*dq54*dq54 + i01*(dq54)*(dq34) + i02*dq34*dq34;
    const double is2 = i10*dq23*dq23 + i11*(dq43)*(dq23) + i12*dq43*dq43;
    const double is3 = i20*dq12*dq12 + i21*(dq32)*(dq12) + i22*dq32*dq32;

    double w1,w2,w3;

    if(WT!=0)
    weno_weights_zt(WT,ct,psi,c1,c2,c3,is1,is2,is3,w1,w2,w3);

    else
    {
    const double is1_psi = is1 + psi;
    const double is2_psi = is2 + psi;
    const double is3_psi = is3 + psi;

    const double a1 = is1_psi*is1_psi;
    const double a2 = is2_psi*is2_psi;
    const double a3 = is3_psi*is3_psi;

    const double sum = c1/a1 + c2/a2 + c3/a3;

    w1 = c1/(epsilon + a1*sum);
    w2 = c2/(epsilon + a2*sum);
    w3 = c3/(epsilon + a3*sum);
    }

    if(SG<0)
    return w1*(q4 + f00*(q3-q4) - f01*(q5-q4))

         + w2*(q3 + f10*(q4-q3) - f11*(q2-q3))

         + w3*(q2 + f20*(q1-q2) + f21*(q3-q2));

    return w1*(q4 + f00*(q3-q4) + f01*(q5-q4))

         + w2*(q3 + f10*(q2-q3) - f11*(q4-q3))

         + w3*(q2 + f20*(q3-q2) - f21*(q1-q2));
}

// One row of x-faces (all j of one i): the stencils run across the rows r1..r6 = f(i-2..i+3,.), the coefficients
// are those of the face, the same for the whole row. No branches: the loop vectorises.
template<int WT>
static inline __attribute__((always_inline)) void sflow_weno_row_x(const int n, const double *__restrict__ r1,
                               const double *__restrict__ r2, const double *__restrict__ r3,
                               const double *__restrict__ r4, const double *__restrict__ r5,
                               const double *__restrict__ r6, double *__restrict__ L, double *__restrict__ R,
                               double ****isf, double ***cf, double ****qf, const int ip,
                               const double psi, const double epsilon, const double ct)
{
    double **const is = isf[ip][1];
    double **const q = qf[ip][1];
    const double *const c = cf[ip][1];

    const double a00=is[0][0], a01=is[0][1], a02=is[0][2], a10=is[1][0], a11=is[1][1], a12=is[1][2],
                 a20=is[2][0], a21=is[2][1], a22=is[2][2];
    const double b00=is[3][0], b01=is[3][1], b02=is[3][2], b10=is[4][0], b11=is[4][1], b12=is[4][2],
                 b20=is[5][0], b21=is[5][1], b22=is[5][2];
    const double c0=c[0], c1=c[1], c2=c[2], c3=c[3], c4=c[4], c5=c[5];
    const double f00=q[0][0], f01=q[0][1], f10=q[1][0], f11=q[1][1], f20=q[2][0], f21=q[2][1];
    const double g00=q[3][0], g01=q[3][1], g10=q[4][0], g11=q[4][1], g20=q[5][0], g21=q[5][1];

    #pragma omp simd
    for(int n1=0; n1<n; ++n1)
    {
    L[n1] = sflow_weno5<WT,-1>(r1[n1],r2[n1],r3[n1],r4[n1],r5[n1], a00,a01,a02,a10,a11,a12,a20,a21,a22,
                               c0,c1,c2, f00,f01,f10,f11,f20,f21, psi,epsilon,ct);

    R[n1] = sflow_weno5<WT,1>(r2[n1],r3[n1],r4[n1],r5[n1],r6[n1], b00,b01,b02,b10,b11,b12,b20,b21,b22,
                              c3,c4,c5, g00,g01,g10,g11,g20,g21, psi,epsilon,ct);
    }
}

// One row of y-faces (all j of one i): the stencil runs along the row r = f(i,.), the coefficients of the face j
// come from the flat tables yt (ytab_fill): coefficient c of face j at yt[c*ny + j].
template<int WT>
static inline __attribute__((always_inline)) void sflow_weno_row_y(const int n, const int ny, const double *__restrict__ r,
                               double *__restrict__ L, double *__restrict__ R, const double *__restrict__ yt,
                               const double psi, const double epsilon, const double ct)
{
    #define YT(c) yt[(c)*ny + n1]
    #pragma omp simd
    for(int n1=0; n1<n; ++n1)
    {
    L[n1] = sflow_weno5<WT,-1>(r[n1-2],r[n1-1],r[n1],r[n1+1],r[n1+2],
                               YT(0),YT(1),YT(2),YT(3),YT(4),YT(5),YT(6),YT(7),YT(8),
                               YT(9),YT(10),YT(11), YT(12),YT(13),YT(14),YT(15),YT(16),YT(17), psi,epsilon,ct);

    R[n1] = sflow_weno5<WT,1>(r[n1-1],r[n1],r[n1+1],r[n1+2],r[n1+3],
                              YT(18),YT(19),YT(20),YT(21),YT(22),YT(23),YT(24),YT(25),YT(26),
                              YT(27),YT(28),YT(29), YT(30),YT(31),YT(32),YT(33),YT(34),YT(35), psi,epsilon,ct);
    }
    #undef YT
}

// The coefficients of the y-faces of the WENO tables (weno_nug_func: isfy, cfy, qfy of vf = 1), 36 per face in
// the order of the arguments of sflow_weno5, left state then right state, one array per coefficient. The tables
// are filled once (weno_nug_func::ini / own_ini) and never change.
void sflow_reconstruct_weno::ytab_fill(lexer *p)
{
    const size_t ny = size_t(MAX(p->knoy,0));
    ytab.assign(36*ny,0.0);

    for(size_t n=0; n<ny; ++n)
    {
    const int jp = int(n) + marge;

        for(int s=0; s<2; ++s)
        {
            for(int a=0; a<3; ++a)
            for(int c=0; c<3; ++c)
            ytab[size_t(18*s + 3*a + c)*ny + n] = isfy[jp][1][3*s+a][c];

            for(int a=0; a<3; ++a)
            ytab[size_t(18*s + 9 + a)*ny + n] = cfy[jp][1][3*s+a];

            for(int a=0; a<3; ++a)
            for(int c=0; c<2; ++c)
            ytab[size_t(18*s + 12 + 2*a + c)*ny + n] = qfy[jp][1][3*s+a][c];
        }
    }
}

// Left and right states of the x-faces of the SLICELOOP1 cells. The WENO values of a whole row of faces are
// computed without branches (sflow_weno_row_x), then stored where the face is active; next to dry cells
// (prepare) the cell value or the central stencil replaces them. The same values as the cell by cell
// reconstruction with weno_nug_func::is_min_x / weight_min_x.
void sflow_reconstruct_weno::reconstruct_x(lexer* p, ghostcell *pgc, fdm2D*, slice& f, slice &fs, slice &fn)
{
    uf=1;
    vf=0;
    wf=0;

    const int nj = p->jmax;
    const int ny = p->knoy;
    const int ie = p->knox - p->ulast;

    if(ny<=0 || ie<=0)
    return;

    if(int(rowL.size())<ny)
    {
    rowL.resize(ny);
    rowR.resize(ny);
    }

    double *const L = rowL.data();
    double *const R = rowR.data();
    const double ps = psi, ep = epsilon, ct = teno_ct;

    for(int ii=0; ii<ie; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;   // IJ of (ii,0)
    const double *const r3 = f.V + q0;
    const double *const r1 = r3 - 2*nj, *const r2 = r3 - nj, *const r4 = r3 + nj, *const r5 = r3 + 2*nj, *const r6 = r3 + 3*nj;
    const int ip = ii + marge;

    if(wtype==1)
    sflow_weno_row_x<1>(ny,r1,r2,r3,r4,r5,r6,L,R,isfx,cfx,qfx,ip,ps,ep,ct);

    else
    if(wtype==2)
    sflow_weno_row_x<2>(ny,r1,r2,r3,r4,r5,r6,L,R,isfx,cfx,qfx,ip,ps,ep,ct);

    else
    sflow_weno_row_x<0>(ny,r1,r2,r3,r4,r5,r6,L,R,isfx,cfx,qfx,ip,ps,ep,ct);

    const int *const flag = p->flagslice1 + q0;
    double *const ws = fs.V + q0;
    double *const wn = fn.V + q0;

        if(!drs)
        {
        for(int jj=0; jj<ny; ++jj)
        if(flag[jj]>0)
        {
        ws[jj] = L[jj];
        wn[jj] = R[jj];
        }
        }

        else
        {
        // dry cells in the stencils of the two cells of the face (prepare)
        const unsigned char *const dl = drx.data() + q0;
        const unsigned char *const dr = dl + nj;
        double **const qf = qfx[ip][1];

        for(int jj=0; jj<ny; ++jj)
        if(flag[jj]>0)
        {
            if(dl[jj]==1)
            ws[jj] = r3[jj];

            else
            if(dl[jj]==2)
            ws[jj] = r3[jj] + qf[1][0]*(r4[jj]-r3[jj]) - qf[1][1]*(r2[jj]-r3[jj]);

            else
            ws[jj] = L[jj];

            if(dr[jj]==1)
            wn[jj] = r4[jj];

            else
            if(dr[jj]==2)
            wn[jj] = r4[jj] + qf[4][0]*(r3[jj]-r4[jj]) - qf[4][1]*(r5[jj]-r4[jj]);

            else
            wn[jj] = R[jj];
        }
        }
    }
}

// Left and right states of the y-faces of the SLICELOOP2 cells, as reconstruct_x.
void sflow_reconstruct_weno::reconstruct_y(lexer* p, ghostcell *pgc, fdm2D*, slice& f, slice &fe, slice &fw)
{
    uf=0;
    vf=1;
    wf=0;

    if(p->j_dir!=1)
    return;

    const int nj = p->jmax;
    const int ny = p->knoy;
    const int je = ny - p->vlast;

    if(je<=0 || p->knox<=0)
    return;

    if(int(ytab.size())!=36*p->knoy)
    ytab_fill(p);

    if(int(rowL.size())<je)
    {
    rowL.resize(je);
    rowR.resize(je);
    }

    double *const L = rowL.data();
    double *const R = rowR.data();
    const double *const yt = ytab.data();
    const double ps = psi, ep = epsilon, ct = teno_ct;

    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;   // IJ of (ii,0)
    const double *const r = f.V + q0;

    if(wtype==1)
    sflow_weno_row_y<1>(je,ny,r,L,R,yt,ps,ep,ct);

    else
    if(wtype==2)
    sflow_weno_row_y<2>(je,ny,r,L,R,yt,ps,ep,ct);

    else
    sflow_weno_row_y<0>(je,ny,r,L,R,yt,ps,ep,ct);

    const int *const flag = p->flagslice2 + q0;
    double *const we = fe.V + q0;
    double *const ww = fw.V + q0;

        if(!drs)
        {
        for(int jj=0; jj<je; ++jj)
        if(flag[jj]>0)
        {
        we[jj] = L[jj];
        ww[jj] = R[jj];
        }
        }

        else
        {
        const unsigned char *const dl = dry.data() + q0;

        for(int jj=0; jj<je; ++jj)
        if(flag[jj]>0)
        {
        double **const qf = qfy[jj+marge][1];

            if(dl[jj]==1)
            we[jj] = r[jj];

            else
            if(dl[jj]==2)
            we[jj] = r[jj] + qf[1][0]*(r[jj+1]-r[jj]) - qf[1][1]*(r[jj-1]-r[jj]);

            else
            we[jj] = L[jj];

            if(dl[jj+1]==1)
            ww[jj] = r[jj+1];

            else
            if(dl[jj+1]==2)
            ww[jj] = r[jj+1] + qf[4][0]*(r[jj]-r[jj+1]) - qf[4][1]*(r[jj+2]-r[jj+1]);

            else
            ww[jj] = R[jj];
        }
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
    const int nj = p->jmax;
    const double *const depth = b->depth.V;
    const double a244 = p->A244;
    const int *const flag1 = p->flagslice1, *const flag2 = p->flagslice2;
    const int ie = p->knox - p->ulast;
    const int je = p->knoy - p->vlast;
    double *const dfx = b->dfx.V, *const dfy = b->dfy.V;
    
    // face depths below the still water level
    for(int ii=0; ii<ie; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag1[q]>0)
        dfx[q] = 0.5*(depth[q+nj]+depth[q]);
    }
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+je; ++q)
        if(flag2[q]>0)
        dfy[q] = 0.5*(depth[q+1]+depth[q]);
    }
    
    pgc->gcsl_start1(p,b->dfx,1);
    pgc->gcsl_start2(p,b->dfy,1);
    
    // face water depths
    {
    const double *const es = b->ETAs.V, *const en = b->ETAn.V;
    double *const Ds = b->Ds.V, *const Dn = b->Dn.V;
    
    for(int ii=0; ii<ie; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag1[q]>0)
        {
        Ds[q] = MAX(es[q] + dfx[q], a244);
        Dn[q] = MAX(en[q] + dfx[q], a244);
        }
    }
    }
    
    {
    const double *const ee = b->ETAe.V, *const ew = b->ETAw.V;
    double *const De = b->De.V, *const Dw = b->Dw.V;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+je; ++q)
        if(flag2[q]>0)
        {
        De[q] = MAX(ee[q] + dfy[q], a244);
        Dw[q] = MAX(ew[q] + dfy[q], a244);
        }
    }
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

