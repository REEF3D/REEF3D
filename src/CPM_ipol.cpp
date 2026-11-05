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

#include"CPM.h"
#include"lexer.h"
#include<limits>

/*--------------------------------------------------------------------
interpolation of the cell-centred fields at the parcels with stencil reuse (WP9)

The parcel update interpolates a dozen cell-centred fields (pressure and stress gradients, solid
fraction, solid velocity, eddy diffusivity, ...) at the same point, and the friction alternates
between the parcel and its substrate point. ccipol4a finds the cell and the weights anew for
every field. cip4a keeps the stencil (cell index and weights) of the last two points and
evaluates the field with it; the arithmetic is that of ccipol4a (lint_a, lint_a_2D), so the
results are identical.
--------------------------------------------------------------------*/

double CPM::cip4a(lexer *p, field &f, double xp, double yp, double zp)
{
    int q;
    
    if(st[0].ok && xp==st[0].x && yp==st[0].y && zp==st[0].z)
    q=0;
    
    else if(st[1].ok && xp==st[1].x && yp==st[1].y && zp==st[1].z)
    q=1;
    
    else
    {
        q = st_last = 1-st_last;
        stencil &s = st[q];
        
        s.x = xp;
        s.y = yp;
        s.z = zp;
        
        s.i = pfi(p,n,xp);
        s.j = pfj(p,n,yp);
        s.k = pfk(p,n,zp);
        
        const int ii=s.i, jj=s.j, kk=s.k;
        
        s.wa = (p->XP[ii+1+marge]-xp)/p->DXP[ii+marge];
        s.wb = (p->YP[jj+1+marge]-yp)/p->DYP[jj+marge];
        s.wc = (p->ZP[kk+1+marge]-zp)/p->DZP[kk+marge];
        s.ok = true;
    }
    
    const stencil &s = st[q];
    const int ii=s.i, kk=s.k;
    const double wa=s.wa, wb=s.wb, wc=s.wc;
    
    if(p->j_dir==0)
    {
        double x1 = wa*f(ii,0,kk)   + (1.0-wa)*f(ii+1,0,kk);
        double x2 = wa*f(ii,0,kk+1) + (1.0-wa)*f(ii+1,0,kk+1);
        
        return wc*x1 +(1.0-wc)*x2;
    }
    
    const int jj=s.j;
    
    double x1 = wa*f(ii,jj,kk)   + (1.0-wa)*f(ii+1,jj,kk);
    double x2 = wa*f(ii,jj+1,kk) + (1.0-wa)*f(ii+1,jj+1,kk);

    double x3 = wa*f(ii,jj,kk+1)   + (1.0-wa)*f(ii+1,jj,kk+1);
    double x4 = wa*f(ii,jj+1,kk+1) + (1.0-wa)*f(ii+1,jj+1,kk+1);

    double y1 = wb*x1 +(1.0-wb)*x2;
    double y2 = wb*x3 +(1.0-wb)*x4;

    return wc*y1 + (1.0-wc)*y2;
}

/*--------------------------------------------------------------------
cell indices of the parcel coordinates (item 4 of the performance review, validation 22)

index_cache stores posc_* and posf_* of the parcel coordinates with the coordinate itself; it runs
at the start of every sub-step and recomputes only the coordinates that changed (in a resting bed
almost none). pci ... pfk return the stored index when the coordinate is the stored one, otherwise
they compute it. The index is a function of the coordinate alone, so the results are those of
posc_* and posf_* whatever the parcel slot q holds (exchanged, removed or added parcels included).
--------------------------------------------------------------------*/
void CPM::index_cache(lexer *p)
{
    const size_t N = P.index;
    
    if(ixX.size()<N)
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        ixX.resize(N,nan);  ixY.resize(N,nan);  ixZ.resize(N,nan);
        ixCi.resize(N,0);   ixCj.resize(N,0);   ixCk.resize(N,0);
        ixFi.resize(N,0);   ixFj.resize(N,0);   ixFk.resize(N,0);
    }
    
    for(size_t q=0;q<N;++q)
    if(P.Flag[q]>=ACTIVE)
    {
        if(P.X[q]!=ixX[q])
        {
            ixX[q] = P.X[q];
            ixCi[q] = p->posc_i(P.X[q]);
            ixFi[q] = p->posf_i(P.X[q]);
        }
        
        if(P.Y[q]!=ixY[q])
        {
            ixY[q] = P.Y[q];
            ixCj[q] = p->posc_j(P.Y[q]);
            ixFj[q] = p->posf_j(P.Y[q]);
        }
        
        if(P.Z[q]!=ixZ[q])
        {
            ixZ[q] = P.Z[q];
            ixCk[q] = p->posc_k(P.Z[q]);
            ixFk[q] = p->posf_k(P.Z[q]);
        }
    }
}

int CPM::pci(lexer *p, int q, double x)
{
    return (q>=0 && size_t(q)<ixX.size() && x==ixX[q]) ? ixCi[q] : p->posc_i(x);
}

int CPM::pcj(lexer *p, int q, double y)
{
    return (q>=0 && size_t(q)<ixY.size() && y==ixY[q]) ? ixCj[q] : p->posc_j(y);
}

int CPM::pck(lexer *p, int q, double z)
{
    return (q>=0 && size_t(q)<ixZ.size() && z==ixZ[q]) ? ixCk[q] : p->posc_k(z);
}

int CPM::pfi(lexer *p, int q, double x)
{
    return (q>=0 && size_t(q)<ixX.size() && x==ixX[q]) ? ixFi[q] : p->posf_i(x);
}

int CPM::pfj(lexer *p, int q, double y)
{
    return (q>=0 && size_t(q)<ixY.size() && y==ixY[q]) ? ixFj[q] : p->posf_j(y);
}

int CPM::pfk(lexer *p, int q, double z)
{
    return (q>=0 && size_t(q)<ixZ.size() && z==ixZ[q]) ? ixFk[q] : p->posf_k(z);
}
