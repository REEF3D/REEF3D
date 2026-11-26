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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_wall.h"
#include<cmath>

// ---------------------------------------------------------------------------------------------------------------------
// NHFLOW semi-implicit step (A 510 11/12/13): wall friction at every wall face (A 519 2), on one Cartesian grid.
// The walls of a cell are those of the explicit NHFLOW and of the turbulence wall functions (nhflow_wall.h): side
// walls in x and y (a solid neighbour or a closed domain side), the bed, a solid above. The log-law drag of a wall
// face on a tangential velocity component u (nhflow_bcmom::wall_friction),
//
//     d(u)/dt = -n_w |u_t| u / (u+^2 dn),   u+ = 1/kappa max(1, ln(30 z0/ks)),   z0 = max(dn/2, ks/30),   ks = B 57,
//
// is implicit on the diagonal of the column operators, with |u_t| of t^n: x walls on v and w, y walls on u and w, a
// solid above (and a solid below a cell above the bottom layer) on u and v. The bed friction of the bottom layer is
// that of A 519 1 (friction: gamma). The rate of a layer face is the mean of the rates of its two cells.
// The layer drag of the next column operator (ld_face, ld_cell: vcD) is the sum of this and of the VRANS drag.
// U and V of a cell have different walls: the column of V takes its own diagonal and its own r = A^-1 1 (wf_vcol).
// ---------------------------------------------------------------------------------------------------------------------

void nhflow_momentum_SI::wf_ini(lexer *p)
{
    wf = (p->A519==2);
    ld = wf || vr;

    wld.assign(p->knoz,0.0);

    if(!wf)
    return;

    wfU.assign(fN,0.0);
    wfV.assign(fN,0.0);
    wfW.assign(fN,0.0);
}

static inline double wf_rate(double ut, double dn, double ks)
{
    const double kappa=0.4;
    double z0 = 0.5*dn;

    if(ks<=0.0)
    ks=0.0001;

    if(30.0*z0<ks)
    z0=ks/30.0;

    const double uplus = (1.0/kappa)*MAX(1.0,std::log(30.0*(z0/ks)));

    return dn>1.0e-12 ? ut/(uplus*uplus*dn) : 0.0;
}

void nhflow_momentum_SI::wf_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // dt times the drag rates of the cells with the velocities of t^n
    const double dt = p->dt;
    int w[6];

    std::fill(wfU.begin(),wfU.end(),0.0);
    std::fill(wfV.begin(),wfV.end(),0.0);
    std::fill(wfW.begin(),wfW.end(),0.0);

    LOOP
    if(p->DF[IJK]>0 && p->wet[IJ]==1 && d->WL(i,j)>p->A544)
    {
    nhflow_wall_faces(p,i,j,k,w);

    const double U = d->U[IJK], V = d->V[IJK]*p->y_dir, W = d->W[IJK];
    const double dz = p->DZN[KP]*d->WL(i,j);

        // x walls: v, w
        int n = w[0]+w[1];
        if(n>0)
        {
        const double r = dt*n*wf_rate(sqrt(V*V + W*W),p->DXN[IP],p->B57);
        wfV[IJK] += r;
        wfW[IJK] += r;
        }

        // y walls: u, w
        n = w[2]+w[3];
        if(n>0)
        {
        const double r = dt*n*wf_rate(sqrt(U*U + W*W),p->DYN[JP],p->B57);
        wfU[IJK] += r;
        wfW[IJK] += r;
        }

        // a solid below a cell above the bottom layer (the bed of the bottom layer: friction), a solid above: u, v
        n = ((w[4] && k>0) ? 1 : 0) + (w[5] ? 1 : 0);
        if(n>0)
        {
        const double r = dt*n*wf_rate(sqrt(U*U + V*V),dz,p->B57);
        wfU[IJK] += r;
        wfV[IJK] += r;
        }
    }

    exchangeC(p,pgc,wfU.data());
    exchangeC(p,pgc,wfV.data());
}

void nhflow_momentum_SI::ld_face(lexer *p, int dir)
{
    // the column of the layer face (i,j) of direction dir for the next col_coef: the inertia factor (VRANS), the
    // layer drag of VRANS and of the walls
    const size_t q0 = (size_t)(lso + i*lsi + j*lsj);

    vcC = vr ? (dir==1 ? vrFCx.data() : vrFCy.data()) + q0 : nullptr;

    if(!wf)
    {
    vcD = (dir==1 ? vrFDx.data() : vrFDy.data()) + q0;
    return;
    }

    const double *W = dir==1 ? wfU.data() : wfV.data();
    const double *V = vr ? (dir==1 ? vrFDx.data() : vrFDy.data()) + q0 : nullptr;
    const int kk = k;

    for(k=0;k<p->knoz;++k)
    wld[k] = (V!=nullptr ? V[k] : 0.0) + 0.5*(W[IJK] + (dir==1 ? W[Ip1JK] : W[IJp1K]));

    k = kk;
    vcD = wld.data();
}

void nhflow_momentum_SI::ld_cell(lexer *p, int comp)
{
    // the column of cell (i,j) for the next col_coef, comp 0 UH, 1 VH, 2 WH
    const size_t q0 = (size_t)(lso + i*lsi + j*lsj);

    vcC = vr ? vrCc.data() + q0 : nullptr;

    if(!wf)
    {
    vcD = (comp==2 ? vrDw.data() : vrDh.data()) + q0;
    return;
    }

    const double *W = comp==0 ? wfU.data() : (comp==1 ? wfV.data() : wfW.data());
    const double *V = vr ? (comp==2 ? vrDw.data() : vrDh.data()) + q0 : nullptr;
    const int kk = k;

    for(k=0;k<p->knoz;++k)
    wld[k] = (V!=nullptr ? V[k] : 0.0) + W[IJK];

    k = kk;
    vcD = wld.data();
}

void nhflow_momentum_SI::wf_vcol(lexer *p, const double *cu, const double *cl, double *cd, double dtgam, double *r)
{
    // the column of VH of cell (i,j) from that of UH (cell_update): the walls of v on the diagonal instead of those of
    // u, and r = A^-1 1 (VRANS: A^-1 C) of that column
    const int kk = k;

    for(k=0;k<p->knoz;++k)
    cd[k] += wfV[IJK] - wfU[IJK];

    k = kk;

    for(int q=0;q<p->knoz;++q)
    wv1[q] = vcC!=nullptr ? vcC[q] : 1.0;

    col_solve(p,cu,cl,cd,dtgam,wv1.data(),r);
}
