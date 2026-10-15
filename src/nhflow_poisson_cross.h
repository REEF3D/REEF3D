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

#ifndef NHFLOW_POISSON_CROSS_H_
#define NHFLOW_POISSON_CROSS_H_

// Sigma cross terms 2 sigx d2P/dx dsig + 2 sigy d2P/dy dsig of the NHFLOW Poisson operator,
// shared by nhflow_poisson and nhflow_poisson_pcorr.
//
// A 529 0 (default): explicit, the cross terms go into the right-hand side with P of the last solve.
//                    This lagged iteration diverges at steep free surfaces (sigx ~ 10, e.g. a run-up
//                    front on a breakwater slope): the NH pressure grows to 100-400 kPa.
// A 529 1:           implicit, the couplings to (i+-1,k+-1), (j+-1,k+-1) are assembled into M.nt..M.eb
//                    wherever the neighbour is an unknown of the solve. REEFMG (N 10 1) multiplies with the
//                    full 15-point (2D: 9-point) operator in its BiCGStab, preconditioned by the V-cycle on
//                    the 7-point part (as FNPF A 328 1). Couplings to non-unknowns (free surface, walls,
//                    dry or shallow columns, below the bed) stay in the right-hand side with the current P.
//                    Only for N 10 1 on the Cartesian grid without mesh refinement, multi-block grids or
//                    membranes; otherwise A 529 1 falls back to 0.
//
// A 530 0 (default): no cross terms in porous cells (POR < 1: VRANS structures, porous bodies, DEM).
// A 530 1:           cross terms also in porous cells.

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"
#include<algorithm>

namespace nhflow_poisson_cross
{

// A 529 1 usable for this assembly
inline bool implicit(lexer *p, fdm_nhf *d)
{
    return p->A529==1 && p->N10==1 && p->geo_curv==0 && p->G1==0 && p->geo_mb==0 && d->MBETA==nullptr;
}

// A 530: cross terms in this row (pressure node between the cells k-1 and k)
inline bool active(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    const int marge = increment::marge;
    
    if(p->A530==1)
    return true;

    return !(d->POR[IJK]<1.0 || d->POR[IJKm1]<1.0);
}

// size and zero the cross couplings for the implicit assembly
inline void ini(lexer *p, fdm_nhf *d, bool impl)
{
    if(!impl)
    {
        // a solve without the cross couplings (fall back): zero couplings, REEFMG's 15-point
        // product then equals the 7-point one
        if(!d->M.nt.empty())
        {
        std::fill(d->M.nt.begin(),d->M.nt.end(),0.0); std::fill(d->M.nb.begin(),d->M.nb.end(),0.0);
        std::fill(d->M.st.begin(),d->M.st.end(),0.0); std::fill(d->M.sb.begin(),d->M.sb.end(),0.0);
        std::fill(d->M.wt.begin(),d->M.wt.end(),0.0); std::fill(d->M.wb.begin(),d->M.wb.end(),0.0);
        std::fill(d->M.et.begin(),d->M.et.end(),0.0); std::fill(d->M.eb.begin(),d->M.eb.end(),0.0);
        }
    return;
    }

    const size_t nr = d->M.p.size();

    d->M.nt.assign(nr,0.0); d->M.nb.assign(nr,0.0);
    d->M.st.assign(nr,0.0); d->M.sb.assign(nr,0.0);
    d->M.wt.assign(nr,0.0); d->M.wb.assign(nr,0.0);
    d->M.et.assign(nr,0.0); d->M.eb.assign(nr,0.0);
}

// pressure node (a,b,kk) is an unknown of the solve: wet, deep, fluid and below the free surface
inline bool unknown(lexer *p, int a, int b, int kk)
{
    if(kk<0 || kk>p->knoz-1)
    return false;

    const int fq = (a-p->imin)*p->jmax*p->kmaxF + (b-p->jmin)*p->kmaxF + kk-p->kmin;
    const int q  = (a-p->imin)*p->jmax + (b-p->jmin);

    return p->flag7[fq]>0 && p->wet[q]==1 && p->deep[q]==1;
}

// implicit: cross couplings of row n with coefficient cX, non-unknown neighbours into the rhs
inline void assemble(lexer *p, fdm_nhf *d, double *P, int i, int j, int k, int n, double cX)
{
    const int marge = increment::marge;
    
    const double ax = cX*2.0*p->sigx[FIJK]/(p->W1*(p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]));

    // rhs form: + ax*(P(i+1,k+1) - P(i-1,k+1) - P(i+1,k-1) + P(i-1,k-1)), matrix = - rhs coefficient
    if(unknown(p,i+1,j,k+1)) d->M.nt[n] = -ax;  else d->rhsvec.V[n] += ax*P[FIp1JKp1];
    if(unknown(p,i-1,j,k+1)) d->M.st[n] =  ax;  else d->rhsvec.V[n] -= ax*P[FIm1JKp1];
    if(unknown(p,i+1,j,k-1)) d->M.nb[n] =  ax;  else d->rhsvec.V[n] -= ax*P[FIp1JKm1];
    if(unknown(p,i-1,j,k-1)) d->M.sb[n] = -ax;  else d->rhsvec.V[n] += ax*P[FIm1JKm1];

    if(p->j_dir==1)
    {
    const double ay = cX*2.0*p->sigy[FIJK]/(p->W1*(p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]));

    if(unknown(p,i,j+1,k+1)) d->M.wt[n] = -ay;  else d->rhsvec.V[n] += ay*P[FIJp1Kp1];
    if(unknown(p,i,j-1,k+1)) d->M.et[n] =  ay;  else d->rhsvec.V[n] -= ay*P[FIJm1Kp1];
    if(unknown(p,i,j+1,k-1)) d->M.wb[n] =  ay;  else d->rhsvec.V[n] -= ay*P[FIJp1Km1];
    if(unknown(p,i,j-1,k-1)) d->M.eb[n] = -ay;  else d->rhsvec.V[n] += ay*P[FIJm1Km1];
    }
}

}

#endif
