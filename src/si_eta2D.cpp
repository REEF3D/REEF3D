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

#include"si_eta2D.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include"solver2D.h"

si_eta2D::si_eta2D(lexer *p) : iter(0), time(0.0), M(p), xvec(p), rhsvec(p)
{
}

si_eta2D::~si_eta2D()
{
}

void si_eta2D::solve(lexer *p, ghostcell *pgc, solver2D *psolv, slice &ax, slice &ay, slice &rhs, slice &eta)
{
    // reaction 1 everywhere: the linear eta system of the fully wet domain
    n=0;
    SLICELOOP4
    {
    assemble(p,ax,ay,rhs,eta,1.0);
    ++n;
    }
    
    run(p,pgc,psolv,eta);
}

void si_eta2D::solve(lexer *p, ghostcell *pgc, solver2D *psolv, slice &ax, slice &ay, slice &rhs, slice &eta, slice &r)
{
    n=0;
    SLICELOOP4
    {
    assemble(p,ax,ay,rhs,eta,r(i,j));
    ++n;
    }
    
    run(p,pgc,psolv,eta);
}

void si_eta2D::assemble(lexer *p, slice &ax, slice &ay, slice &rhs, slice &eta, double r)
{
    // row n of cell (i,j), SLICELOOP4 order
    double cn,cs,cw,ce;
    
    cn = ax(i,j)/p->DXN[IP];
    cs = ax(i-1,j)/p->DXN[IP];
    cw = 0.0;
    ce = 0.0;
    
    if(p->j_dir==1)
    {
    cw = ay(i,j)/p->DYN[JP];
    ce = ay(i,j-1)/p->DYN[JP];
    }
    
    M.p[n] = r + cn + cs + cw + ce;
    M.n[n] = -cn;
    M.s[n] = -cs;
    M.w[n] = -cw;
    M.e[n] = -ce;
    
    rhsvec.V[n] = rhs(i,j);
    
    // Dirichlet: a face towards a cell outside the domain (open water level boundary)
    if(cn!=0.0 && p->flagslice4[Ip1J]<0)
    {
    rhsvec.V[n] += cn*eta(i+1,j);
    M.n[n] = 0.0;
    }
    
    if(cs!=0.0 && p->flagslice4[Im1J]<0)
    {
    rhsvec.V[n] += cs*eta(i-1,j);
    M.s[n] = 0.0;
    }
    
    if(cw!=0.0 && p->flagslice4[IJp1]<0)
    {
    rhsvec.V[n] += cw*eta(i,j+1);
    M.w[n] = 0.0;
    }
    
    if(ce!=0.0 && p->flagslice4[IJm1]<0)
    {
    rhsvec.V[n] += ce*eta(i,j-1);
    M.e[n] = 0.0;
    }
    
    // no coupling and no reaction (cannot happen for r = 1): keep eta
    if(M.p[n]<=0.0)
    {
    M.p[n] = 1.0;
    rhsvec.V[n] = eta(i,j);
    }
}

void si_eta2D::run(lexer *p, ghostcell *pgc, solver2D *psolv, slice &eta)
{
    const double starttime = pgc->timer();
    
    psolv->start(p,pgc,eta,M,xvec,rhsvec,4);
    
    iter = p->solveriter;
    time = pgc->timer()-starttime;
}
