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
#include<vector>
#include<algorithm>
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include"solver2D.h"

si_eta2D::si_eta2D(lexer *p) : iter(0), time(0.0), M(p), xvec(p), rhsvec(p), pw(p), pw0(p), pr(p), rW(p), react(p), rhs(p)
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
    {
    pw(i,j) = 1.0;
    pw0(i,j) = 1.0;
    pr(i,j) = 0.0;
    }
}

void si_eta2D::wd_rows(lexer *p, slice &rhsW, slice &depth, slice &bdiag)
{
    SLICELOOP4
    {
    react(i,j) = pw(i,j) + bdiag(i,j);
    rhs(i,j) = rhsW(i,j) - pw(i,j)*depth(i,j);
    rW(i,j) = rhsW(i,j);
    }
}

int si_eta2D::wd_update(lexer *p, slice &ax, slice &ay, slice &bdiag, slice &eta, slice &depth, bool islands)
{
    int change=0;
    
    SLICELOOP4
    {
        const bool coupled = ax(i,j)>0.0 || ax(i-1,j)>0.0 || bdiag(i,j)>0.0
                          || (p->j_dir==1 && (ay(i,j)>0.0 || ay(i,j-1)>0.0));
        
        double pn = (eta(i,j) + depth(i,j)>0.0 || !coupled) ? 1.0 : 0.0;
        
        pw(i,j) = pn;
    }
    
    // dry islands with water become wet (wd_islands)
    wd_nisl = wd_nisl_cells = wd_nisl_empty = 0;
    
    if(islands)
    wd_islands(p,ax,ay,bdiag);
    
    SLICELOOP4
    {
        if(pw(i,j)!=pw0(i,j))
        ++change;
        
    pw0(i,j) = pw(i,j);
    }
    
    return change;
}

int si_eta2D::wd_islands(lexer *p, slice &ax, slice &ay, slice &bdiag)
{
    // A set of dry rows (pw = 0, eta + d <= 0) whose open faces lead only to each other has no
    // reaction, no boundary and no exchange with the rest of the domain: its rows are a pure-Neumann
    // block of the system (any constant solves the homogeneous part, the volume balance
    // sum(rhs) = 0 is not met in general), the solution drifts there and the multigrid coarse
    // operator of the block is singular (zero coarse diagonal, reefmg_core::coarsen_t). Such an island
    // has nowhere to drain to, so it keeps its water: if it holds water (the sum of its right-hand sides,
    // its volume, is positive) its cells become wet again (pw = 1) and the Newton iteration finds the
    // level that holds that volume; an empty island stays dry. The search is local to the grid of this rank: a face
    // towards a cell outside it (another rank, another block, a boundary) counts as an exit, so only
    // islands that lie wholly inside the grid of the rank are found.
    const int nx = p->knox, ny = p->knoy;
    
    if((int)wdmark.size()!=nx*ny)
    wdmark.assign(nx*ny,0);
    else
    std::fill(wdmark.begin(),wdmark.end(),0);
    
    auto inside = [&](int a, int b) { return a>=0 && b>=0 && a<nx && b<ny; };
    
    // the open neighbours of (a,b): x faces ax(a,b) (i+1), ax(a-1,b) (i-1), y faces ay(a,b) (j+1), ay(a,b-1) (j-1)
    const int di[4] = {1,-1,0,0}, dj[4] = {0,0,1,-1};
    auto open = [&](int a, int b, int d) -> bool
    {
        if(d==0) return ax(a,b)>0.0;
        if(d==1) return ax(a-1,b)>0.0;
        if(p->j_dir!=1) return false;
        if(d==2) return ay(a,b)>0.0;
        return ay(a,b-1)>0.0;
    };
    
    int nisl=0, ncell=0, nempty=0;
    std::vector<int> comp, stack;
    
    SLICELOOP4
    {
        const int q0 = i*ny + j;
        
        if(pw(i,j)>0.5 || wdmark[q0]!=0)
        continue;
        
        // the component of the dry rows connected to (i,j) through open faces
        comp.clear();
        stack.clear();
        stack.push_back(q0);
        wdmark[q0]=1;
        bool exit_ = false;
        
        while(!stack.empty())
        {
            const int q = stack.back();
            stack.pop_back();
            comp.push_back(q);
            const int a = q/ny, b = q%ny;
            
            if(bdiag(a,b)>0.0)
            exit_ = true;
            
            for(int d=0; d<4; ++d)
            {
                if(!open(a,b,d))
                continue;
                
                const int a2 = a+di[d], b2 = b+dj[d];
                
                if(!inside(a2,b2))
                {
                exit_ = true;
                continue;
                }
                
                const int i0=i, j0=j;
                i=a2; j=b2;
                const bool cell = p->flagslice4[IJ]>0;
                i=i0; j=j0;
                
                if(!cell || pw(a2,b2)>0.5)
                {
                exit_ = true;
                continue;
                }
                
                const int q2 = a2*ny + b2;
                
                if(wdmark[q2]==0)
                {
                wdmark[q2]=1;
                stack.push_back(q2);
                }
            }
        }
        
        if(!exit_)
        {
            // the water of the island: the sum of its volume equations (the face terms cancel), and
            // whether it was made wet before in this step
            double vol = 0.0, area = 0.0;
            bool again = false;
            for(int q : comp)
            {
            i = q/ny;
            j = q%ny;
            const double A = p->geo_curv==1 ? p->AREA2D[IJ] : p->DXN[IP]*(p->j_dir==1 ? p->DYN[JP] : 1.0);
            vol += rW(i,j)*A;
            area += A;
            again = again || pr(i,j)>0.5;
            }
            i = q0/ny;
            j = q0%ny;

            // An island with water has wet cells: all wet as the start of the Newton iteration on the
            // island, the next passes dry the cells above its level again (the cells below it keep the
            // island connected to wet cells). An island without water to speak of (mean depth below
            // wd_film: the solve does not resolve its level, its cells would dry again) stays dry, as
            // does one that was made wet before in this step (at most once per step: no ping-pong);
            // its equations are (nearly) consistent and the multigrid keeps a singular coarse block
            // finite (reefmg_core::coarsen_t).
            if(vol<=wd_film*area || again)
            {
            ++nempty;
            continue;
            }

            for(int q : comp)
            {
            pw(q/ny,q%ny) = 1.0;
            pr(q/ny,q%ny) = 1.0;
            }
            
        ++nisl;
        ncell += (int)comp.size();
        }
    }
    
    wd_nisl = nisl;
    wd_nisl_cells = ncell;
    wd_nisl_empty = nempty;
    
    return nisl;
}
