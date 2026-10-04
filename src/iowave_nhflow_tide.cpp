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

void iowave::nhflow_open_edges(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH, slice &WL)
{
    if(p->open_xm==0 && p->open_xp==0)
    return;
    
    if(edge_h.empty())
    edge_h.assign(size_t(p->imax)*size_t(p->jmax),0.0);
    
    const double g = fabs(p->W22);
    
    for(int side : {1,2})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr)
        continue;
        
        const int b = bgs.index(z->bg);
        const int count = side==1 ? p->gcin_count : p->gcout_count;
        int **gc = side==1 ? p->gcin : p->gcout;
        for(n=0;n<count;++n)
        {
        i=gc[n][0];
        j=gc[n][1];
        k=gc[n][2];
        
        if(gc[n][3]!=(side==1 ? 1 : 4))
        continue;
        
            const double h0 = d->depth(i,j);
            const double hi = d->WL(i,j);
            
            if(p->wet[IJ]==0 || hi<=1.0e-6 || h0<=1.0e-6)
            continue;
            
            const double eb = bgs.eta(b);
            double ub, vb;
            bgs.vel(b,h0,ub,vb);
            
            const double ui = nhflow_col_ubar(p,d,d->U);
            double ug, hg;
            
            if(z->method==bc_method::riemann)
            {
                const double hb = fmax(h0+eb,1.0e-6);
                double Rin, Rout;
                
                if(side==1)
                {
                Rin  = ub + 2.0*sqrt(g*hb);
                Rout = ui - 2.0*sqrt(g*hi);
                hg = pow(Rin-Rout,2.0)/(16.0*g);
                }
                else
                {
                Rin  = ub - 2.0*sqrt(g*hb);
                Rout = ui + 2.0*sqrt(g*hi);
                hg = pow(Rout-Rin,2.0)/(16.0*g);
                }
                
                ug = 0.5*(Rin+Rout);
            }
            else
            {
                // outward normal: -x at x-, +x at x+
                hg = hi;
                ug = side==1 ? ub - sqrt(g/hi)*((hi-h0) - eb) : ub + sqrt(g/hi)*((hi-h0) - eb);
            }
            
            edge_h[IJ] = hg;
            
            // inflow: tangential velocity of the background, outflow: of the interior
            const bool in = (side==1) ? ug>0.0 : ug<0.0;
            
            const double vg = in ? vb : V[IJK];
            const double wg = in ? 0.0 : W[IJK];
            
            if(side==1)
            {
            U[Im1JK]=U[Im2JK]=U[Im3JK]=ug;
            V[Im1JK]=V[Im2JK]=V[Im3JK]=vg;
            W[Im1JK]=W[Im2JK]=W[Im3JK]=wg;
            UH[Im1JK]=UH[Im2JK]=UH[Im3JK]=hg*ug;
            VH[Im1JK]=VH[Im2JK]=VH[Im3JK]=hg*vg;
            WH[Im1JK]=WH[Im2JK]=WH[Im3JK]=hg*wg;
            }
            else
            {
            U[Ip1JK]=U[Ip2JK]=U[Ip3JK]=ug;
            V[Ip1JK]=V[Ip2JK]=V[Ip3JK]=vg;
            W[Ip1JK]=W[Ip2JK]=W[Ip3JK]=wg;
            UH[Ip1JK]=UH[Ip2JK]=UH[Ip3JK]=hg*ug;
            VH[Ip1JK]=VH[Ip2JK]=VH[Ip3JK]=hg*vg;
            WH[Ip1JK]=WH[Ip2JK]=WH[Ip3JK]=hg*wg;
            }
        }
    }
}

void iowave::nhflow_open_edges_rk(lexer *p, fdm_nhf *d, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    // x+ ghost cells of the RK arrays (the x- ones are copied for all inflow cells)
    if(p->open_xp==0)
    return;
    
    for(n=0;n<p->gcout_count;++n)
    {
    i=p->gcout[n][0];
    j=p->gcout[n][1];
    k=p->gcout[n][2];
    
    if(p->gcout[n][3]!=4)
    continue;
    
        U[Ip1JK]=d->U[Ip1JK]; U[Ip2JK]=d->U[Ip2JK]; U[Ip3JK]=d->U[Ip3JK];
        V[Ip1JK]=d->V[Ip1JK]; V[Ip2JK]=d->V[Ip2JK]; V[Ip3JK]=d->V[Ip3JK];
        W[Ip1JK]=d->W[Ip1JK]; W[Ip2JK]=d->W[Ip2JK]; W[Ip3JK]=d->W[Ip3JK];
        UH[Ip1JK]=d->UH[Ip1JK]; UH[Ip2JK]=d->UH[Ip2JK]; UH[Ip3JK]=d->UH[Ip3JK];
        VH[Ip1JK]=d->VH[Ip1JK]; VH[Ip2JK]=d->VH[Ip2JK]; VH[Ip3JK]=d->VH[Ip3JK];
        WH[Ip1JK]=d->WH[Ip1JK]; WH[Ip2JK]=d->WH[Ip2JK]; WH[Ip3JK]=d->WH[Ip3JK];
    }
}

void iowave::nhflow_open_edges_wl(lexer *p, fdm_nhf *d, slice &WL)
{
    if(p->open_xm==0 && p->open_xp==0)
    return;
    
    for(int side : {1,2})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr)
        continue;
        
        const int count = side==1 ? p->gcslin_count : p->gcslout_count;
        int **gc = side==1 ? p->gcslin : p->gcslout;
        
        for(n=0;n<count;++n)
        {
        i=gc[n][0];
        j=gc[n][1];
        
        if(gc[n][3]!=(side==1 ? 1 : 4))
        continue;
        
            // Riemann: h_g of this step; Flather: zero gradient
            const double hg = (z->method==bc_method::riemann && !edge_h.empty() && edge_h[IJ]>0.0) ? edge_h[IJ] : WL(i,j);
            const int s = side==1 ? -1 : 1;
            
            for(int q=1; q<=3; ++q)
            {
            WL(i+s*q,j) = hg;
            d->eta(i+s*q,j) = hg - d->depth(i,j);
            }
        }
    }
}
