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


#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"

/*--------------------------------------------------------------------
periodic boundaries for the parcels (DIVEMesh C 21 x, C 22 y)

serial (periodic = 1): the kernel deposits into the ghost cells beyond the periodic sides,
pfold adds them to the cells on the opposite side before start4a_sum (which copies the
periodic ghost cells and would overwrite the deposits). Parcels crossing a periodic side are
moved by the period after the grid-limited step, so the limiter sees the move into the ghost
cell, which is the copy of the cell on the opposite side.

parallel (periodic = 2): the ranks at the periodic sides are neighbours, start4a_sum sums the
deposits and the parcels go to the neighbour by the regular exchange; the receiving rank moves
them by the period (part::xchange_fillback_flag).
--------------------------------------------------------------------*/

void CPM::pfold(lexer *p, field &f)
{
    if(perx==1)
    {
        for(j=-1;j<p->knoy+1;++j)
        for(k=-1;k<p->knoz+1;++k)
        {
            f(p->knox-1,j,k) += f(-1,j,k);
            f(0,j,k) += f(p->knox,j,k);
            f(-1,j,k) = 0.0;
            f(p->knox,j,k) = 0.0;
        }
    }
    
    if(pery==1)
    {
        for(i=-1;i<p->knox+1;++i)
        for(k=-1;k<p->knoz+1;++k)
        {
            f(i,p->knoy-1,k) += f(i,-1,k);
            f(i,0,k) += f(i,p->knoy,k);
            f(i,-1,k) = 0.0;
            f(i,p->knoy,k) = 0.0;
        }
    }
}

// serial periodic: parcels whose new position PX1 lies beyond a periodic side are moved by the
// period, together with their old position PX0 (the RK2 stages combine both)
void CPM::periodic_wrap(lexer *p, double *PX0, double *PY0, double *PX1, double *PY1)
{
    double Lx = p->global_xmax - p->global_xmin;
    double Ly = p->global_ymax - p->global_ymin;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        if(perx==1)
        {
            if(PX1[n]>=p->global_xmax)
            {
                PX1[n] -= Lx;
                PX0[n] -= Lx;
            }
            else if(PX1[n]<p->global_xmin)
            {
                PX1[n] += Lx;
                PX0[n] += Lx;
            }
        }
        
        if(pery==1)
        {
            if(PY1[n]>=p->global_ymax)
            {
                PY1[n] -= Ly;
                PY0[n] -= Ly;
            }
            else if(PY1[n]<p->global_ymin)
            {
                PY1[n] += Ly;
                PY0[n] += Ly;
            }
        }
    }
}

// serial periodic: the ghost cells beyond a periodic side carry the flags of the cells on the
// opposite side, so the interpolation of the fluid velocities (which skips solid faces) uses the
// periodic ghost values next to the periodic side
void CPM::periodic_flags(lexer *p)
{
    auto id = [&](int ii, int jj, int kk) {return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;};
    int *fl[4] = {p->flag1,p->flag2,p->flag3,p->flag4};
    
    if(perx==1)
    for(int jj=p->jmin;jj<p->jmin+p->jmax;++jj)
    for(int kk=p->kmin;kk<p->kmin+p->kmax;++kk)
    for(int q=1;q<=-p->imin;++q)
    for(int f=0;f<4;++f)
    {
        fl[f][id(-q,jj,kk)] = fl[f][id(p->knox-q,jj,kk)];
        fl[f][id(p->knox-1+q,jj,kk)] = fl[f][id(q-1,jj,kk)];
    }
    
    if(pery==1)
    for(int ii=p->imin;ii<p->imin+p->imax;++ii)
    for(int kk=p->kmin;kk<p->kmin+p->kmax;++kk)
    for(int q=1;q<=-p->imin;++q)
    for(int f=0;f<4;++f)
    {
        fl[f][id(ii,-q,kk)] = fl[f][id(ii,p->knoy-q,kk)];
        fl[f][id(ii,p->knoy-1+q,kk)] = fl[f][id(ii,q-1,kk)];
    }
}
