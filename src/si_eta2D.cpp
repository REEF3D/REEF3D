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

si_eta2D::si_eta2D(lexer *p) : iter(0), time(0.0), M(p), xvec(p), rhsvec(p), pw(p), react(p), rhs(p)
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
    
    // curvilinear grid (G 50, NHFLOW): the face lengths over the cell area
    if(p->geo_curv==1)
    {
    const double A = p->AREA2D[IJ];
    
    cn = ax(i,j)*p->LEN1[IJ]/A;
    cs = ax(i-1,j)*p->LEN1[Im1J]/A;
    cw = ay(i,j)*p->LEN2[IJ]/A;
    ce = ay(i,j-1)*p->LEN2[IJm1]/A;
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

int si_eta2D::solve_wetdry(lexer *p, ghostcell *pgc, solver2D *psolv, slice &ax, slice &ay, slice &rhsW, slice &depth, slice &bdiag, slice &eta, int itmax, int &dry)
{
    int it_done=0;
    
    wd_start(p);
    
    for(int it=0; it<itmax; ++it)
    {
    wd_rows(p,rhsW,depth,bdiag);
    
    solve(p,pgc,psolv,ax,ay,rhs,eta,react);
    
        if(p->mpi_size>1)
        {
        pgc->gcslparax(p,eta,4);
        pgc->gcslparacox(p,eta,14);
        }
    
    ++it_done;
    
        if(itmax==1)
        break;
        
        // new wet pattern
        int change = wd_update(p,ax,ay,bdiag,eta,depth);
        
        change = pgc->globalisum(change);
        
        if(change==0)
        break;
    }
    
    dry=0;
    
    SLICELOOP4
    if(pw(i,j)<0.5)
    ++dry;
    
    dry = pgc->globalisum(dry);
    
    return it_done;
}

void si_eta2D::build(lexer *p, slice &ax, slice &ay, slice &rhs, slice &eta, slice &r)
{
    // the rows of the system (M, rhsvec) without solving it (composite solve of a multi-block grid)
    n=0;
    SLICELOOP4
    {
    assemble(p,ax,ay,rhs,eta,r(i,j));
    ++n;
    }
}

void si_eta2D::wd_start(lexer *p)
{
    SLICELOOP4
    pw(i,j) = 1.0;
}

void si_eta2D::wd_rows(lexer *p, slice &rhsW, slice &depth, slice &bdiag)
{
    SLICELOOP4
    {
    react(i,j) = pw(i,j) + bdiag(i,j);
    rhs(i,j) = rhsW(i,j) - pw(i,j)*depth(i,j);
    }
}

int si_eta2D::wd_update(lexer *p, slice &ax, slice &ay, slice &bdiag, slice &eta, slice &depth)
{
    int change=0;
    
    SLICELOOP4
    {
        const bool coupled = ax(i,j)>0.0 || ax(i-1,j)>0.0 || bdiag(i,j)>0.0
                          || (p->j_dir==1 && (ay(i,j)>0.0 || ay(i,j-1)>0.0));
        
        const double pn = (eta(i,j) + depth(i,j)>0.0 || !coupled) ? 1.0 : 0.0;
        
        if(pn!=pw(i,j))
        ++change;
        
    pw(i,j) = pn;
    }
    
    return change;
}
