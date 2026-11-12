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

#include"si_tridiag2D.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include"matrix2D.h"
#include"vec2D.h"

si_tridiag2D::si_tridiag2D(lexer *p, ghostcell *pgc, solver2D *fallback) : fb(fallback)
{
}

si_tridiag2D::~si_tridiag2D()
{
    delete fb;
}

void si_tridiag2D::start(lexer *p, ghostcell *pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var)
{
    // the rows in SLICELOOP4 order and their i; row n couples to n-1 (M.s) and n+1 (M.n) only, the neighbours in i
    ci.clear();
    SLICELOOP4
    ci.push_back(i);

    const int N = (int)ci.size();
    bool tri = true;

    for(int n=0; n<N && tri; ++n)
    {
        if(M.w[n]!=0.0 || M.e[n]!=0.0)
        tri = false;

        if(M.n[n]!=0.0 && (n+1>=N || ci[n+1]!=ci[n]+1))
        tri = false;

        if(M.s[n]!=0.0 && (n==0 || ci[n-1]!=ci[n]-1))
        tri = false;
    }

    if(!tri)
    {
    fb->start(p,pgc,f,M,xvec,rhsvec,var);
    return;
    }

    // Thomas: a_n x_{n-1} + b_n x_n + c_n x_{n+1} = r_n (diagonally dominant: b = r + cn + cs >= |a| + |c|)
    c.resize(N);
    d.resize(N);

    for(int n=0; n<N; ++n)
    {
        const double a = n>0 ? M.s[n] : 0.0;
        const double den = M.p[n] - (n>0 ? a*c[n-1] : 0.0);
        c[n] = M.n[n]/den;
        d[n] = (rhsvec.V[n] - (n>0 ? a*d[n-1] : 0.0))/den;
    }

    for(int n=N-2; n>=0; --n)
    d[n] -= c[n]*d[n+1];

    int n=0;
    SLICELOOP4
    {
    f(i,j) = d[n];
    ++n;
    }

    p->solveriter = 1;
    p->final_res = 0.0;
}
