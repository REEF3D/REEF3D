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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

/*--------------------------------------------------------------------
Tidal / current background in NHFLOW (iowave redesign, step 4):

- relaxation zones with a background (B 523) target background + waves,
  the wave ramp acts on the waves only; beach zones with a background relax
  to the background instead of still water (iowave_nhflow_relax.cpp)
- Riemann edge (B 520 method 3): ghost cells from the incoming characteristic
  of the background and the outgoing one of the interior,
    x-: R+ = u_b + 2 sqrt(g h_b), R- = u_i - 2 sqrt(g h_i)
    u_g = (R+ + R-)/2, h_g = (R+ - R-)^2/(16 g)
- Flather edge (B 520 method 4): h_g = h_i, u_g = u_b -+ sqrt(g/h_i) (eta_i - eta_b)
  (x- / x+), i.e. q_n = q_b + sqrt(g h) (eta - eta_b) with the outward normal

The edges set U, V, W, UH, VH, WH (depth uniform) and WL, eta in the three
ghost cells; ghostcell treats the edge as open (lexer open_xm / open_xp).
--------------------------------------------------------------------*/

void iowave::bg_build(lexer *p)
{
    col_gen_bg.assign(size_t(p->imax)*size_t(p->jmax),-1);
    col_beach_bg.assign(size_t(p->imax)*size_t(p->jmax),-1);
    
    ILOOP
    JLOOP
    {
        const double x = p->XP[IP];
        const double y = p->YP[JP];
        
        const bc_zone *zg = zones.relax_zone_at(x,y);
        if(zg!=nullptr && zg->bg>0)
        col_gen_bg[IJ] = bgs.index(zg->bg);
        
        const bc_zone *zb = zones.beach_zone_at(x,y);
        if(zb!=nullptr && zb->bg>0)
        col_beach_bg[IJ] = bgs.index(zb->bg);
    }
    
    bg_built = true;
}

int iowave::gen_bg(lexer *p)
{
    return col_gen_bg[IJ];
}

int iowave::beach_bg(lexer *p)
{
    return col_beach_bg[IJ];
}

void iowave::nhflow_bg_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!bg_on)
    return;
    
    if(!bg_built)
    bg_build(p);
    
    bgs.update(p,p->simtime);
    
    // still water depth per column
    col_h0.assign(size_t(p->imax)*size_t(p->jmax),0.0);
    
    SLICELOOP4
    col_h0[IJ] = d->depth(i,j);
}

double iowave::nhflow_col_ubar(lexer *p, fdm_nhf *d, double *F)
{
    // depth average of a layer field in the column (i,j); sigma spacing DZN sums to 1
    const int kk = k;
    double s = 0.0;
    
    for(k=0; k<p->knoz; ++k)
    s += F[IJK]*p->DZN[KP];
    
    k = kk;
    return s;
}

// side code of the ghostcell lists (1: x-, 2: y+, 3: y-, 4: x+) and ghost offset of an edge (1: x-, 2: x+, 3: y-, 4: y+)
static void edge_geometry(int e, int &sc, int &di, int &dj, double &nx, double &ny)
{
    sc = e==1 ? 1 : e==2 ? 4 : e==3 ? 3 : 2;
    di = e==1 ? -1 : e==2 ? 1 : 0;
    dj = e==3 ? -1 : e==4 ? 1 : 0;
    nx = -double(di);    // inward normal
    ny = -double(dj);
}

// edge e is a Riemann / Flather edge with its ghost cells set here (lexer flags, set in iowave::ini)
static int edge_open(const lexer *p, int e)
{
    return e==1 ? p->open_xm : e==2 ? p->open_xp : e==3 ? p->open_ym : p->open_yp;
}

void iowave::nhflow_open_edges(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH, slice &WL)
{
    if(p->open_xm==0 && p->open_xp==0 && p->open_ym==0 && p->open_yp==0)
    return;
    
    if(edge_h.empty())
    edge_h.assign(size_t(p->imax)*size_t(p->jmax),0.0);
    
    const double g = fabs(p->W22);
    
    for(int side : {1,2,3,4})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr || edge_open(p,side)==0)
        continue;
        
        int sc, di, dj;
        double nx, ny;
        edge_geometry(side,sc,di,dj,nx,ny);
        const bool xedge = side<=2;
        
        const int b = bgs.index(z->bg);
        
        // waves of the edge's sources (B 524, Riemann only): eta and depth averaged u, v per column
        const bool waves = z->method==bc_method::riemann && !z->sources.empty();
        
        if(waves)
        {
            select_sources(&z->sources);
            
            if(edge_etaw.empty())
            {
            edge_etaw.assign(size_t(p->imax)*size_t(p->jmax),0.0);
            edge_uw.assign(size_t(p->imax)*size_t(p->jmax),0.0);
            edge_vw.assign(size_t(p->imax)*size_t(p->jmax),0.0);
            }
            
            for(int list=0; list<2; ++list)
            {
            const int cs = list==0 ? p->gcslin_count : p->gcslout_count;
            int **gs = list==0 ? p->gcslin : p->gcslout;
            
            for(n=0;n<cs;++n)
            {
            i=gs[n][0];
            j=gs[n][1];
            
            if(gs[n][3]!=sc)
            continue;
            
                xg = xgen(p);
                yg = ygen(p);
                
                edge_etaw[IJ] = wave_eta(p,pgc,xg,yg);
                
                double su = 0.0, sv = 0.0;
                for(k=0; k<p->knoz; ++k)
                {
                su += wave_u(p,pgc,xg,yg,p->ZSP[IJK]-p->phimean)*p->DZN[KP];
                
                if(!xedge)
                sv += wave_v(p,pgc,xg,yg,p->ZSP[IJK]-p->phimean)*p->DZN[KP];
                }
                
                edge_uw[IJ] = su;
                edge_vw[IJ] = sv;
            }
            }
        }
        
        for(int list=0; list<2; ++list)
        {
        const int count = list==0 ? p->gcin_count : p->gcout_count;
        int **gc = list==0 ? p->gcin : p->gcout;
        
        for(n=0;n<count;++n)
        {
        i=gc[n][0];
        j=gc[n][1];
        k=gc[n][2];
        
        if(gc[n][3]!=sc)
        continue;
        
            const double h0 = d->depth(i,j);
            const double hi = d->WL(i,j);
            
            if(p->wet[IJ]==0 || hi<=1.0e-6 || h0<=1.0e-6)
            continue;
            
            // background at the boundary face
            const double xf = side==1 ? p->XN[IP] : side==2 ? p->XN[IP1] : p->XP[IP];
            const double yf = side==3 ? p->YN[JP] : side==4 ? p->YN[JP1] : p->YP[JP];
            const double eb = bgs.eta(b,xf,yf);
            double ub, vb;
            bgs.vel(b,h0,xf,yf,ub,vb);
            
            // normal (inward) and tangential velocities: x edges U / V, y edges V / U
            const double nn = xedge ? nx : ny;
            const double ui = xedge ? nhflow_col_ubar(p,d,d->U) : nhflow_col_ubar(p,d,d->V);
            
            // waves of the layer and of the column (ramped), added to the background
            double uwk=0.0, vwk=0.0, wwk=0.0, ewc=0.0, uwc=0.0, vwc=0.0;
            
            if(waves)
            {
                xg = xgen(p);
                yg = ygen(p);
                const double zk = p->ZSP[IJK]-p->phimean;
                const double rw = ramp(p);
                
                uwk = rw*wave_u(p,pgc,xg,yg,zk);
                vwk = rw*wave_v(p,pgc,xg,yg,zk);
                wwk = rw*wave_w(p,pgc,xg,yg,zk);
                ewc = rw*edge_etaw[IJ];
                uwc = rw*edge_uw[IJ];
                vwc = rw*edge_vw[IJ];
            }
            
            const double unb = xedge ? nn*(ub + uwc) : nn*(vb + vwc);   // target, normal
            const double uni = nn*ui;                                  // interior, normal
            double ung, hg;
            
            if(z->method==bc_method::riemann)
            {
                const double hb = fmax(h0+eb+ewc,1.0e-6);
                const double Rin  = unb + 2.0*sqrt(g*hb);
                const double Rout = uni - 2.0*sqrt(g*hi);
                
                hg  = pow(Rin-Rout,2.0)/(16.0*g);
                ung = 0.5*(Rin+Rout);
            }
            else
            {
                // q_n(outward) = q_n,b + sqrt(g h) (eta - eta_b)
                hg = hi;
                ung = (xedge ? nn*ub : nn*vb) - sqrt(g/hi)*((hi-h0) - eb);
            }
            
            edge_h[IJ] = hg;
            
            // inflow: tangential velocity of the background (+ waves), outflow: of the interior
            const bool in = ung>0.0;
            
            // depth average from the characteristics, vertical profile of the waves
            double ug, vg;
            
            if(xedge)
            {
            ug = nn*ung + (uwk - uwc);
            vg = in ? vb + vwk : V[IJK];
            }
            else
            {
            vg = nn*ung + (vwk - vwc);
            ug = in ? ub + uwk : U[IJK];
            }
            
            const double wg = in ? wwk : W[IJK];
            
            const int ii=i, jj=j;
            for(int q=1; q<=3; ++q)
            {
                i = ii + q*di;
                j = jj + q*dj;
                
                U[IJK] = ug;
                V[IJK] = vg;
                W[IJK] = wg;
                UH[IJK] = hg*ug;
                VH[IJK] = hg*vg;
                WH[IJK] = hg*wg;
            }
            i = ii;
            j = jj;
        }
        }
        
        if(waves)
        select_sources(nullptr);
    }
}

void iowave::nhflow_open_edges_rk(lexer *p, fdm_nhf *d, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    // ghost cells of the RK arrays from the step's ghost cells
    for(int side : {1,2,3,4})
    {
        if(edge_open(p,side)==0)
        continue;
        
        int sc, di, dj;
        double nx, ny;
        edge_geometry(side,sc,di,dj,nx,ny);
        
        for(int list=0; list<2; ++list)
        {
        const int count = list==0 ? p->gcin_count : p->gcout_count;
        int **gc = list==0 ? p->gcin : p->gcout;
        
        for(n=0;n<count;++n)
        {
        i=gc[n][0];
        j=gc[n][1];
        k=gc[n][2];
        
        if(gc[n][3]!=sc)
        continue;
        
            const int ii=i, jj=j;
            for(int q=1; q<=3; ++q)
            {
                i = ii + q*di;
                j = jj + q*dj;
                
                U[IJK]=d->U[IJK];
                V[IJK]=d->V[IJK];
                W[IJK]=d->W[IJK];
                UH[IJK]=d->UH[IJK];
                VH[IJK]=d->VH[IJK];
                WH[IJK]=d->WH[IJK];
            }
            i = ii;
            j = jj;
        }
        }
    }
}

void iowave::nhflow_open_edges_wl(lexer *p, fdm_nhf *d, slice &WL)
{
    if(p->open_xm==0 && p->open_xp==0 && p->open_ym==0 && p->open_yp==0)
    return;
    
    for(int side : {1,2,3,4})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr || edge_open(p,side)==0)
        continue;
        
        int sc, di, dj;
        double nx, ny;
        edge_geometry(side,sc,di,dj,nx,ny);
        
        for(int list=0; list<2; ++list)
        {
        const int count = list==0 ? p->gcslin_count : p->gcslout_count;
        int **gc = list==0 ? p->gcslin : p->gcslout;
        
        for(n=0;n<count;++n)
        {
        i=gc[n][0];
        j=gc[n][1];
        
        if(gc[n][3]!=sc)
        continue;
        
            // Riemann: h_g of this step; Flather: zero gradient
            const double hg = (z->method==bc_method::riemann && !edge_h.empty() && edge_h[IJ]>0.0) ? edge_h[IJ] : WL(i,j);
            
            for(int q=1; q<=3; ++q)
            {
            WL(i+q*di,j+q*dj) = hg;
            d->eta(i+q*di,j+q*dj) = hg - d->depth(i,j);
            }
        }
        }
    }
}
