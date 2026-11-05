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

#include"sflow_amr_linesolve.h"
#include"lexer.h"
#include"slice.h"
#include"matrix2D.h"
#include"vec2D.h"
#include<cmath>

void sflow_amr_linesolve::start(lexer *p, ghostcell *pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhs, int var)
{
    const int nsl = p->imax*p->jmax;
    auto L = [&](int ii, int jj) { return (ii-p->imin)*p->jmax + (jj-p->jmin); };

    id.assign(nsl,-1);
    bool xline=true, yline=true;

    n=0;
    SLICELOOP4
    {
    id[L(i,j)] = n;
    if(M.w[n]!=0.0 || M.e[n]!=0.0) xline=false;
    if(M.n[n]!=0.0 || M.s[n]!=0.0) yline=false;
    ++n;
    }

    // one line: cells k0..k1 along the direction, Thomas algorithm
    auto line = [&](int len, auto cell)
    {
        int k=0;
        while(k<len)
        {
            int ii,jj;
            cell(k,ii,jj);
            if(id[L(ii,jj)]<0) { ++k; continue; }

            int k0=k;
            while(k<len) { cell(k,ii,jj); if(id[L(ii,jj)]<0) break; ++k; }
            int m = k-k0;

            a.resize(m); bb.resize(m); c.resize(m); d.resize(m);
            for(int q=0; q<m; ++q)
            {
                cell(k0+q,ii,jj);
                int r = id[L(ii,jj)];
                a[q]  = xline ? M.s[r] : M.e[r];
                c[q]  = xline ? M.n[r] : M.w[r];
                bb[q] = M.p[r];
                d[q]  = rhs.V[r];
            }

            // ends: neighbours outside the computed cells are Dirichlet values
            {
                int i0,j0,i1,j1;
                cell(k0,i0,j0);
                cell(k0+m-1,i1,j1);
                if(xline) { d[0] -= a[0]*f(i0-1,j0); d[m-1] -= c[m-1]*f(i1+1,j1); }
                else      { d[0] -= a[0]*f(i0,j0-1); d[m-1] -= c[m-1]*f(i1,j1+1); }
            }

            for(int q=1; q<m; ++q)
            {
                double w = a[q]/bb[q-1];
                bb[q] -= w*c[q-1];
                d[q]  -= w*d[q-1];
            }
            d[m-1] /= bb[m-1];
            for(int q=m-2; q>=0; --q)
            d[q] = (d[q] - c[q]*d[q+1])/bb[q];

            for(int q=0; q<m; ++q)
            {
                cell(k0+q,ii,jj);
                f(ii,jj) = d[q];
            }
        }
    };

    if(xline)
    {
        for(int jj=0; jj<p->knoy; ++jj)
        line(p->knox,[&](int k, int &ii, int &j2){ ii=k; j2=jj; });
    }
    else if(yline)
    {
        for(int ii=0; ii<p->knox; ++ii)
        line(p->knoy,[&](int k, int &i2, int &jj){ i2=ii; jj=k; });
    }
    else
    {
        // general 5-point rows: symmetric Gauss-Seidel
        for(int it=0; it<200; ++it)
        {
            double dmax=0.0;
            for(int pass=0; pass<2; ++pass)
            for(int s=0; s<p->knox; ++s)
            for(int t=0; t<p->knoy; ++t)
            {
                int ii = pass==0 ? s : p->knox-1-s;
                int jj = pass==0 ? t : p->knoy-1-t;
                int r = id[L(ii,jj)];
                if(r<0) continue;
                double v = (rhs.V[r] - M.n[r]*f(ii+1,jj) - M.s[r]*f(ii-1,jj) - M.w[r]*f(ii,jj+1) - M.e[r]*f(ii,jj-1))/M.p[r];
                dmax = fmax(dmax,fabs(v-f(ii,jj)));
                f(ii,jj) = v;
            }
            if(dmax<1.0e-12)
            break;
        }
    }

    p->solveriter = 1;
}
