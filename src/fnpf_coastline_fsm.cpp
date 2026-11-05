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


#include"fnpf_coastline.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

// Coastline signed distance by fast sweeping with closest-point propagation
// (8-neighbour vector distance transform, Danielsson 1980; sweep orders as in Zhao 2005).
//
// Replaces the pseudo-time PDE reinitialisation at startup, which needs
// 2*maxlength/DXM RK3 iterations (3 WENO passes + 3 halo swaps each).
// Here every cell carries the coordinates of its closest interface point, and
// alternating Gauss-Seidel sweeps pass these on to the neighbours. Halo swap
// after each sweep, repeated until no value changes on any rank.
// The result is the Euclidean distance to the interface points, independent of
// the domain decomposition. Works on non-uniform grids.
//
// Input : f = +1 (wet) / -1 (dry), as set in fnpf_coastline::start
// Output: f = signed distance to the wet/dry interface; the interface is placed
//         at the cell-face midpoint between a wet and a dry cell.
//
//    L    : unsigned distance
//    cpx  : x-coordinate of the closest interface point
//    cpy  : y-coordinate of the closest interface point
//    frk1 : sign (+1 wet, -1 dry)
//    frk2 : 1 = frozen (interface or inactive cell), 0 = to be swept

void fnpf_coastline::fsm(lexer *p, ghostcell *pgc, slice &f)
{
    const double big = 1.0e20;
    const double dmax = 2.0*p->maxlength;   // cap for cells no interface reaches (no coastline in domain)
    double dxh,dyh,xs,ys,xc,yc,nx,ny,t,diff,Lold;
    int sweep,iter,maxiter,nsweep;
    
    if(p->mpirank==0)
    cout<<"initializing coastline (fast sweeping)... "<<endl<<endl;
    
    pgc->gcsl_start4(p,f,50);
    
    // ------------------------------------------------------------
    // 1. seeds: interface cells get their closest point on the wet/dry
    //    interface through the adjacent face midpoints, and are frozen
    // ------------------------------------------------------------
    SLICEBASELOOP
    {
    L(i,j) = big;
    cpx(i,j) = big;
    cpy(i,j) = big;
    frk1(i,j) = 1.0;
    frk2(i,j) = 1.0;   // inactive cells: never updated
    }
    
    SLICELOOP4
    {
        frk1(i,j) = f(i,j)>=0.0?1.0:-1.0;
        frk2(i,j) = 0.0;
        
        xc = p->XP[IP];
        yc = p->YP[JP];
        
        dxh = big;
        dyh = big;
        xs = 0.0;
        ys = 0.0;
        
        // signed offsets to the nearest wet/dry face midpoint in x and y
        if(f(i-1,j)*frk1(i,j)<0.0 && 0.5*p->DXP[IM1]<dxh)
        {
        dxh = 0.5*p->DXP[IM1];
        xs = -dxh;
        }
        
        if(f(i+1,j)*frk1(i,j)<0.0 && 0.5*p->DXP[IP]<dxh)
        {
        dxh = 0.5*p->DXP[IP];
        xs = dxh;
        }
        
        if(p->j_dir==1)
        {
        if(f(i,j-1)*frk1(i,j)<0.0 && 0.5*p->DYP[JM1]<dyh)
        {
        dyh = 0.5*p->DYP[JM1];
        ys = -dyh;
        }
        
        if(f(i,j+1)*frk1(i,j)<0.0 && 0.5*p->DYP[JP]<dyh)
        {
        dyh = 0.5*p->DYP[JP];
        ys = dyh;
        }
        }
        
        // corner cell: project the cell centre onto the line through both midpoints
        if(dxh<big && dyh<big)
        {
        nx = 1.0/xs;
        ny = 1.0/ys;
        t  = 1.0/(nx*nx + ny*ny);   // line: nx*dx + ny*dy = 1
        cpx(i,j) = xc + nx*t;
        cpy(i,j) = yc + ny*t;
        L(i,j) = sqrt(t);
        frk2(i,j) = 1.0;
        }
        
        else
        if(dxh<big)
        {
        cpx(i,j) = xc + xs;
        cpy(i,j) = yc;
        L(i,j) = dxh;
        frk2(i,j) = 1.0;
        }
        
        else
        if(dyh<big)
        {
        cpx(i,j) = xc;
        cpy(i,j) = yc + ys;
        L(i,j) = dyh;
        frk2(i,j) = 1.0;
        }
    }
    
    pgc->gcsl_start4(p,L,50);
    pgc->gcsl_start4(p,cpx,50);
    pgc->gcsl_start4(p,cpy,50);
    
    // ------------------------------------------------------------
    // 2. sweeping
    // ------------------------------------------------------------
    nsweep  = p->j_dir==1?4:2;
    maxiter = MAX(p->M10,1) + 10;   // info crosses at least one subdomain per outer iteration
    
    for(iter=0; iter<maxiter; ++iter)
    {
        diff=0.0;
        
        for(sweep=0; sweep<nsweep; ++sweep)
        {
            const int idir = (sweep%2==0)?1:-1;
            const int jdir = (sweep<2)?1:-1;
            
            for(int ii=0; ii<p->knox; ++ii)
            {
            i = idir>0?ii:p->knox-1-ii;
            
                for(int jj=0; jj<p->knoy; ++jj)
                {
                j = jdir>0?jj:p->knoy-1-jj;
                
                if(frk2(i,j)<0.5)
                {
                Lold = L(i,j);
                
                fsm_update(p,i-1,j);
                fsm_update(p,i+1,j);
                
                if(p->j_dir==1)
                {
                fsm_update(p,i,j-1);
                fsm_update(p,i,j+1);
                fsm_update(p,i-1,j-1);
                fsm_update(p,i+1,j-1);
                fsm_update(p,i-1,j+1);
                fsm_update(p,i+1,j+1);
                }
                
                diff = MAX(diff, Lold-L(i,j));
                }
                }
            }
            
            pgc->gcsl_start4(p,L,50);
            pgc->gcsl_start4(p,cpx,50);
            pgc->gcsl_start4(p,cpy,50);
        }
        
        diff = pgc->globalmax(diff);
        
        if(diff < 1.0e-10*p->DXM)
        break;
    }
    
    if(p->mpirank==0)
    cout<<"coastline fast sweeping iterations: "<<MIN(iter+1,maxiter)<<endl<<endl;
    
    if(p->mpirank==0 && iter==maxiter)
    cout<<"WARNING: coastline fast sweeping not converged"<<endl<<endl;
    
    // ------------------------------------------------------------
    // 3. signed distance
    // ------------------------------------------------------------
    SLICELOOP4
    f(i,j) = frk1(i,j)*MIN(L(i,j),dmax);
    
    pgc->gcsl_start4(p,f,50);
}

// take over the closest interface point of neighbour (in,jn) if it is closer
void fnpf_coastline::fsm_update(lexer *p, int in, int jn)
{
    if(cpx(in,jn) >= 1.0e19)
    return;
    
    const double ddx = p->XP[IP] - cpx(in,jn);
    const double ddy = p->YP[JP] - cpy(in,jn);
    const double d = sqrt(ddx*ddx + ddy*ddy);
    
    if(d < L(i,j))
    {
    L(i,j)   = d;
    cpx(i,j) = cpx(in,jn);
    cpy(i,j) = cpy(in,jn);
    }
}
