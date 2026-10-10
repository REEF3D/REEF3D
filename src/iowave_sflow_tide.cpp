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
#include"fdm2D.h"
#include"ghostcell.h"
#include<iomanip>

/*--------------------------------------------------------------------
Tidal / current backgrounds and open edges in SFLOW (as NHFLOW,
iowave_nhflow_tide.cpp), for the depth-averaged cell-centred HLL scheme:

- relaxation zones with a background (B 523) target background + waves,
  beach zones relax to the background (iowave_sflow_relax.cpp)
- Riemann (B 520 method 3), Flather (4), clamped level (5) and clamped
  discharge (6, B 525) edges set U, V, W and eta in the three ghost cells of
  the in- / outflow lists; sflow_momentum_func::vel_bc leaves open sides alone
  and ghostcell (gcsl_epol4) leaves their eta to iowave
- waves of the edge's sources (B 524) at a Riemann edge: eta and the velocities
  averaged over B 160 + 1 levels, as the SFLOW precalc
- waves on the background (B 530) and the edge mass balance (B 529)

The ghost water level of a Riemann / clamped level edge comes from the step's
characteristics (edge_h) and is set again after every water level update
(waterlevel2D -> sflow_open_edges_eta).
--------------------------------------------------------------------*/

// side code of the ghostcell lists (1: x-, 2: y+, 3: y-, 4: x+) and ghost offset of an edge (1: x-, 2: x+, 3: y-, 4: y+)
static void sflow_edge_geometry(int e, int &sc, int &di, int &dj, double &nx, double &ny)
{
    sc = e==1 ? 1 : e==2 ? 4 : e==3 ? 3 : 2;
    di = e==1 ? -1 : e==2 ? 1 : 0;
    dj = e==3 ? -1 : e==4 ? 1 : 0;
    nx = -double(di);    // inward normal
    ny = -double(dj);
}

static int sflow_edge_open(const lexer *p, int e)
{
    return e==1 ? p->open_xm : e==2 ? p->open_xp : e==3 ? p->open_ym : p->open_yp;
}

void iowave::sflow_bg_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(!bg_on)
    return;

    if(!bg_built)
    bg_build(p);

    bgs.update(p,p->simtime);

    // still water depth per column
    col_h0.assign(size_t(p->imax)*size_t(p->jmax),0.0);

    SLICELOOP4
    col_h0[IJ] = b->depth(i,j);

    // waves on the background (B 530): k on h_eff, Doppler
    if(p->B530>0)
    nhflow_wave_background(p,pgc);

    // edge mass balance (B 529)
    if(p->B529>0)
    sflow_mass_balance(p,b,pgc);
}

void iowave::sflow_open_edges(lexer *p, fdm2D *b, ghostcell *pgc, slice &U, slice &V)
{
    if(p->open_xm==0 && p->open_xp==0 && p->open_ym==0 && p->open_yp==0)
    return;

    if(edge_h.empty())
    edge_h.assign(size_t(p->imax)*size_t(p->jmax),0.0);

    const double g = fabs(p->W22);

    for(int side : {1,2,3,4})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr || sflow_edge_open(p,side)==0)
        continue;

        int sc, di, dj;
        double nx, ny;
        sflow_edge_geometry(side,sc,di,dj,nx,ny);
        const bool xedge = side<=2;

        const int bi = z->bg>0 ? bgs.index(z->bg) : -1;

        // clamped discharge edge with B 525: normal velocity Q / (sum of h ds along the edge), ramped
        double uq = 0.0;

        if(z->method==bc_method::clamp_q && z->has_Q)
        {
            double area = 0.0;

            for(int list=0; list<2; ++list)
            {
            const int cs = list==0 ? p->gcslin_count : p->gcslout_count;
            int **gs = list==0 ? p->gcslin : p->gcslout;

            for(n=0;n<cs;++n)
            if(gs[n][3]==sc)
            {
                i=gs[n][0];
                j=gs[n][1];

                if(p->wet[IJ]==1)
                area += b->WL(i,j)*(xedge ? p->DYN[JP] : p->DXN[IP]);
            }
            }

            area = pgc->globalsum(area);

            const double pi = 3.14159265358979323846;
            const double r = (z->Q_tramp>0.0 && p->simtime<z->Q_tramp) ? 0.5*(1.0-cos(pi*fmax(p->simtime,0.0)/z->Q_tramp)) : 1.0;

            uq = area>1.0e-12 ? r*z->Q/area : 0.0;
        }

        // waves of the edge's sources (B 524, Riemann only)
        const bool waves = z->method==bc_method::riemann && !z->sources.empty();

        if(waves)
        select_sources(&z->sources);

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

            const double h0 = b->depth(i,j);
            const double hi = b->eta(i,j) + b->depth(i,j);

            if(p->wet[IJ]==0 || hi<=1.0e-6 || h0<=1.0e-6)
            continue;

            // background at the boundary face
            const double xf = side==1 ? p->XN[IP] : side==2 ? p->XN[IP1] : p->XP[IP];
            const double yf = side==3 ? p->YN[JP] : side==4 ? p->YN[JP1] : p->YP[JP];
            const double eb = bi>=0 ? bgs.eta(bi,xf,yf) : 0.0;
            double ub=0.0, vb=0.0;
            if(bi>=0)
            bgs.vel(bi,h0,xf,yf,ub,vb);

            const double nn = xedge ? nx : ny;
            const double ui = xedge ? U(i,j) : V(i,j);

            // waves (ramped): eta and velocities averaged over B 160 + 1 levels
            double ewc=0.0, uwc=0.0, vwc=0.0, wwc=0.0;

            if(waves)
            {
                xg = xgen(p);
                yg = ygen(p);
                const double rw = ramp(p);
                const double ew = wave_eta(p,pgc,xg,yg);
                const double dz = (eb + ew + p->wd - b->bed(i,j))/double(p->B160);

                double su=0.0, sv=0.0, sw=0.0;
                double zl = -p->wd;
                for(int qn=0; qn<=p->B160; ++qn)
                {
                su += wave_u(p,pgc,xg,yg,zl);
                if(p->j_dir==1)
                sv += wave_v(p,pgc,xg,yg,zl);
                sw += wave_w(p,pgc,xg,yg,zl);
                zl += dz;
                }

                ewc = rw*ew;
                uwc = rw*su/double(p->B160+1);
                vwc = rw*sv/double(p->B160+1);
                wwc = rw*sw/double(p->B160+1);
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
            if(z->method==bc_method::flather)
            {
                hg = hi;
                ung = (xedge ? nn*ub : nn*vb) - sqrt(g/hi)*((hi-h0) - eb);
            }
            else
            if(z->method==bc_method::clamp_level)
            {
                hg = fmax(h0+eb,1.0e-6);
                ung = uni;
            }
            else
            {
                hg = hi;
                ung = z->has_Q ? uq : (xedge ? nn*ub : nn*vb);
            }

            edge_h[IJ] = hg;

            // inflow: tangential velocity of the background (+ waves), outflow: of the interior
            const bool in = ung>0.0;

            double ug, vg;

            if(xedge)
            {
            ug = nn*ung;
            vg = in ? vb + vwc : V(i,j);
            }
            else
            {
            vg = nn*ung;
            ug = in ? ub + uwc : U(i,j);
            }

            const double wg = in ? wwc : b->W(i,j);
            const double eg = (z->method==bc_method::riemann || z->method==bc_method::clamp_level) ? hg - h0 : b->eta(i,j);

            for(int q=1; q<=3; ++q)
            {
                const int ii = i + q*di;
                const int jj = j + q*dj;

                U(ii,jj) = ug;
                V(ii,jj) = vg;
                b->W(ii,jj) = wg;
                b->UA(ii,jj) = ug;
                b->VA(ii,jj) = vg;
                b->eta(ii,jj) = eg;
            }
        }
        }

        if(waves)
        select_sources(nullptr);
    }
}

void iowave::sflow_open_edges_eta(lexer *p, fdm2D *b, slice &eta)
{
    if(p->open_xm==0 && p->open_xp==0 && p->open_ym==0 && p->open_yp==0)
    return;

    for(int side : {1,2,3,4})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr || sflow_edge_open(p,side)==0)
        continue;

        int sc, di, dj;
        double nx, ny;
        sflow_edge_geometry(side,sc,di,dj,nx,ny);

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

            // Riemann, clamped level: h_g of this step; Flather, clamped discharge: zero gradient
            const bool own_h = z->method==bc_method::riemann || z->method==bc_method::clamp_level;
            const double eg = (own_h && !edge_h.empty() && edge_h[IJ]>0.0) ? edge_h[IJ] - b->depth(i,j) : eta(i,j);

            for(int q=1; q<=3; ++q)
            eta(i+q*di,j+q*dj) = eg;
        }
        }
    }
}

void iowave::sflow_mass_balance(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double t = p->simtime;

    if(t==mb_t)
    return;

    // flux into the domain through each open edge [m3/s], from the continuity fluxes FEx / FEy
    double q[5]={0.0,0.0,0.0,0.0,0.0};

    for(int side : {1,2,3,4})
    {
        const bc_zone *z = zones.open_edge(side);
        if(z==nullptr || sflow_edge_open(p,side)==0)
        continue;

        int sc, di, dj;
        double nx, ny;
        sflow_edge_geometry(side,sc,di,dj,nx,ny);

        for(int list=0; list<2; ++list)
        {
        const int cs = list==0 ? p->gcslin_count : p->gcslout_count;
        int **gs = list==0 ? p->gcslin : p->gcslout;

        for(int m=0; m<cs; ++m)
        if(gs[m][3]==sc)
        {
            i = gs[m][0];
            j = gs[m][1];

            double f = 0.0;
            if(side==1) f =  b->FEx(i-1,j);
            if(side==2) f = -b->FEx(i,j);
            if(side==3) f =  b->FEy(i,j-1);
            if(side==4) f = -b->FEy(i,j);

            q[side] += f*(side<=2 ? p->DYN[JP] : p->DXN[IP]);
        }
        }

        q[side] = pgc->globalsum(q[side]);
    }

    double V = 0.0;
    SLICELOOP4
    V += b->WL(i,j)*p->DXN[IP]*p->DYN[JP];
    V = pgc->globalsum(V);

    if(mb_t<0.0)
    {
        if(p->mpirank==0)
        {
            mb_out.open("./REEF3D_SFLOW_Log/REEF3D-SFLOW-iowave-mass-balance.dat");
            mb_out<<"# iowave edge mass balance (B 529 "<<p->B529<<"): flux into the domain [m3/s], mean over the interval"<<std::endl;
            mb_out<<"# residual = dV/dt - sum of the edges: relaxation / beach zones and time discretisation"<<std::endl;
            mb_out<<"# t  V  dV/dt";
            for(int side : {1,2,3,4})
            if(zones.open_edge(side)!=nullptr && sflow_edge_open(p,side)==1)
            mb_out<<"  Q_edge"<<side<<"(zone "<<zones.open_edge(side)->id<<")";
            mb_out<<"  sum  residual"<<std::endl;
        }

        mb_t = mb_t0 = t;
        mb_V0 = V;
        mb_n = 0;
        return;
    }

    const double dt = t - mb_t;
    for(int side : {1,2,3,4})
    mb_int[side] += q[side]*dt;

    mb_t = t;
    ++mb_n;

    if(mb_n<p->B529)
    return;

    const double T = t - mb_t0;

    if(p->mpirank==0 && T>0.0)
    {
        double sum = 0.0;
        mb_out<<std::setprecision(10)<<t<<"  "<<V<<"  "<<(V-mb_V0)/T;
        for(int side : {1,2,3,4})
        if(zones.open_edge(side)!=nullptr && sflow_edge_open(p,side)==1)
        {
            mb_out<<"  "<<mb_int[side]/T;
            sum += mb_int[side]/T;
        }
        mb_out<<"  "<<sum<<"  "<<(V-mb_V0)/T - sum<<std::endl;
    }

    for(int side : {1,2,3,4})
    mb_int[side] = 0.0;

    mb_t0 = t;
    mb_V0 = V;
    mb_n = 0;
}
