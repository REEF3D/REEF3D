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
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include<vector>

/*--------------------------------------------------------------------
sub-grid bedload layer (Q 58 1, S 10 1)

The bedload layer is a few grains thick, far below the grid scale. The parcels of the
layer are not moved by the resolved drag and the grid-limited step, but by a sub-grid
closure for each bed column, while they still carry the sediment mass (bed evolution,
grading, slopes and exchange with the suspension follow from the parcels):

  Shields vector of a grain d in a bed column, with the bed slope (gravity along the bed):
      theta = tau_b/((rho_s-rho_f) g d) - theta_c,d/mu_s grad(z_b)
      tau_b       : bed shear stress, log law at the first fluid cell at least 0.5 h above the bed
                    (distance from the bed level set, ks = S 21 d50)
      theta_c,d   : Q 60 (d/d50)^(-0.8)  (hiding of the small grains, equal mobility tendency)
      mu_s        : static friction Q 36
      the slope term is limited to tan(beta) = 0.9 mu_s: without flow the layer stays at rest,
      slopes beyond the angle of repose avalanche through the packed-bed stress

  grain velocity (Fernandez Luque & van Beek 1976):
      u_b = Q 59 sqrt(R g d) (sqrt|theta| - 0.7 sqrt(theta_c,d)),  along theta

  moving volume per bed area at equilibrium (Bagnold 1956: the moving grains carry the excess
  shear stress with the dynamic friction mu_d = mu_s):
      V_m* = (|theta| - theta_c,d) d/mu_s

  pickup and deposition (Einstein 1950): a grain is picked up from the exposed top layer, travels
  a hop of length L (exponential distribution, mean Q 62 d) with u_b and is deposited on the
  bed; a grain in a column below the threshold is deposited at once. Pickup rate per bed area
      E = V_m* u_b/L     ->  equilibrium V_m = V_m*,  q = V_m* u_b,
  out of equilibrium the transport adapts over the hop length L. With Q 59 = 6.5, mu_s = 0.63
  the equilibrium rate follows Meyer-Peter & Mueller within 8 % for 0.06 < theta < 0.5
  (Wong & Parker 2006: Q 59 = 3.2).

  The exposed grains of a column (exposure_update) share the pickup of the column by their
  share of the bed area, expo V_parcel/(theta_0 d). A parcel holds many grains, a column holds
  a fraction of a parcel at the equilibrium: the pickup is accumulated per column and a parcel
  is picked up when the accumulated volume reaches the parcel volume.

  Q 58 2, pickup into the suspension: where u* > w_s (Bagnold's criterion), parcels of the
  bedload layer are released into the resolved flow 0.5 h above the bed with the van Rijn (1984)
  reference concentration,
      E_s = w_s c_a,   c_a = 0.015 d50 T^1.5/(a D*^0.3),  a = 0.5 h,  T = theta/theta_c - 1
  the resolved suspension (drag, settling, turbulent dispersion Q 52) brings them back to the bed.

The flow moves the bed through the layer only:
  - the exposed grains are fully sheltered (fluid at rest at the grain level, Q 57 off), grains
    resting within one cell above the bed level see no resolved flow and no turbulent dispersion;
    moving grains near the bed (settling, avalanches) see the resolved flow reduced as delta/h
  - the jumps of the layer parcels (pickup, hops, deposition) bypass the grid-limited step and
    are placed into the first cell with free volume (capacity as in the grid-limited step)
  - bed surfaces (CPM_bedchange.cpp): the layer uses a single valued bed level per column from the
    resting parcels of the top cells (no staircase of the iso-surface); the fluid sees this level
    relaxed in time (Q 63, 10 s), so that single parcels moving in and out of the top cells do not
    switch fluid cells on and off; the parcels keep the iso-surface of the solid fraction for
    exposure, the near-bed closure, the pressure and their own forcing inside the bed
--------------------------------------------------------------------*/

// bed column of a position: periodic sides wrap, beyond the local domain the edge column
void CPM::bedload_column(lexer *p, double x, double y, int &ic, int &jc)
{
    if(perx==1)
    {
        double Lx = p->global_xmax - p->global_xmin;
        if(x>=p->global_xmax) x-=Lx;
        if(x<p->global_xmin) x+=Lx;
    }

    if(pery==1)
    {
        double Ly = p->global_ymax - p->global_ymin;
        if(y>=p->global_ymax) y-=Ly;
        if(y<p->global_ymin) y+=Ly;
    }

    ic = MAX(0, MIN(p->knox-1, p->posc_i(x)));
    jc = p->j_dir==1 ? MAX(0, MIN(p->knoy-1, p->posc_j(y))) : 0;
}

// bed shear stress vector and bed slope per column
//   log law of the wall at the centre of the first fluid cell at least 0.5 h above the bed,
//   with the distance to the bed from the topo level set (consistent with the wall function
//   of the fluid; the linear interpolation towards the solid cells below must not be used)
void CPM::bedload_columns(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    const double ks = MAX(p->S21*p->S20, 1.0e-6);

    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        int kb = MAX(0,MIN(p->knoz-1, p->posc_k(s->bedzh(i,j))));
        double ur=0.0, vr=0.0, zr=0.0;
        bool found=false, rest=false;

        for(k=kb;k<p->knoz;++k)
        {
            double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);

            if(a->topo(i,j,k) >= 0.5*h)
            {
                // a solid body above the bed (e.g. a pipeline with a small gap): the fluid cell
                // below it if there is one, else no flow (the body rests on the bed)
                bool body = p->solidread>0 && a->solid(i,j,k)<0.0;
                
                if(body && k>kb && a->topo(i,j,k-1)>0.0 && a->solid(i,j,k-1)>0.0)
                {
                    --k;
                    body = false;
                }
                
                rest = body;
                ur = body ? 0.0 : 0.5*(a->u(i-1,j,k) + a->u(i,j,k));
                vr = (body || p->j_dir==0) ? 0.0 : 0.5*(a->v(i,j-1,k) + a->v(i,j,k));
                zr = a->topo(i,j,k);
                blH(i,j) = h;
                found = true;
                break;
            }
        }

        // no fluid cell above the bed in this subdomain (vertical decomposition): no transport
        if(!found)
        {
            blTx(i,j) = blTy(i,j) = 0.0;
            blH(i,j) = p->DZN[kb+marge];
            continue;
        }

        // reference height Q 66 (in cells above the bed, > 0): the velocity of the log law is taken at
        // Q 66 h above the bed, linear between the cell centres of the column, at most half way to
        // a solid body or the free surface above; the first cell next to the bed is under-resolved
        // in an accelerated flow (gap below a pipeline) and gives a too low bed shear stress
        if(p->Q66>0.0 && !rest)
        {
            int k0 = k;
            double h = blH(i,j);
            double zbed = p->ZP[k0+marge] - a->topo(i,j,k0);
            double zref = zbed + p->Q66*h;
            
            // top of the fluid column above the bed: solid body or free surface
            for(int kk=k0+1; kk<p->knoz; ++kk)
            {
                double dtop = 1.0e20;
                
                if(p->solidread>0 && a->solid(i,j,kk)<0.0)
                dtop = p->ZP[kk+marge] + a->solid(i,j,kk) - zbed;
                
                if(a->phi(i,j,kk)<0.0)
                dtop = MIN(dtop, p->ZP[kk+marge] + a->phi(i,j,kk) - zbed);
                
                if(dtop<1.0e19)
                {
                    zref = MIN(zref, zbed + 0.5*dtop);
                    break;
                }
            }
            
            zref = MAX(zref, zbed + zr);
            
            // linear between the cell centres around zref
            int kr = k0;
            while(kr+1<p->knoz && p->ZP[kr+1+marge]<=zref)
            ++kr;
            
            if(kr>k0 && kr+1<p->knoz && !(p->solidread>0 && a->solid(i,j,kr+1)<0.0))
            {
                double f = (zref - p->ZP[kr+marge])/(p->ZP[kr+1+marge] - p->ZP[kr+marge]);
                double u0 = 0.5*(a->u(i-1,j,kr) + a->u(i,j,kr)), u1 = 0.5*(a->u(i-1,j,kr+1) + a->u(i,j,kr+1));
                ur = (1.0-f)*u0 + f*u1;
                
                if(p->j_dir==1)
                {
                    double v0 = 0.5*(a->v(i,j-1,kr) + a->v(i,j,kr)), v1 = 0.5*(a->v(i,j-1,kr+1) + a->v(i,j,kr+1));
                    vr = (1.0-f)*v0 + f*v1;
                }
                
                zr = zref - zbed;
            }
            else if(kr>k0)
            {
                ur = 0.5*(a->u(i-1,j,kr) + a->u(i,j,kr));
                vr = p->j_dir==1 ? 0.5*(a->v(i,j-1,kr) + a->v(i,j,kr)) : 0.0;
                zr = p->ZP[kr+marge] - zbed;
            }
        }

        double um = sqrt(ur*ur + vr*vr);
        double uplus = log(MAX(30.0*zr/ks, 1.0+1.0e-6))/0.4;
        double us = um/uplus;
        double tau = p->W1*us*us;

        blTx(i,j) = um>1.0e-12 ? tau*ur/um : 0.0;
        blTy(i,j) = um>1.0e-12 ? tau*vr/um : 0.0;
    }

    // bed slope, central differences; periodic sides wrap, walls one-sided
    auto zb = [&](int ii, int jj, bool &ok) -> double
    {
        ok = true;

        if(ii<0 || ii>=p->knox)
        {
            if(perx==1)
            ii = ii<0 ? ii+p->knox : ii-p->knox;

            else if((ii<0 && p->nb1<0) || (ii>=p->knox && p->nb4<0))
            ok = false;
        }

        if(jj<0 || jj>=p->knoy)
        {
            if(pery==1)
            jj = jj<0 ? jj+p->knoy : jj-p->knoy;

            else if((jj<0 && p->nb3<0) || (jj>=p->knoy && p->nb2<0))
            ok = false;
        }

        return ok ? zbl(s,ii,jj) : 0.0;
    };

    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        bool okm,okp;
        double zc = zbl(s,i,j);
        double zm = zb(i-1,j,okm);
        double zp = zb(i+1,j,okp);
        double dx = (okm ? p->DXP[IM1] : 0.0) + (okp ? p->DXP[IP] : 0.0);

        blGx(i,j) = dx>0.0 ? ((okp ? zp : zc) - (okm ? zm : zc))/dx : 0.0;

        blGy(i,j) = 0.0;

        if(p->j_dir==1)
        {
            zm = zb(i,j-1,okm);
            zp = zb(i,j+1,okp);
            double dy = (okm ? p->DYP[JM1] : 0.0) + (okp ? p->DYP[JP] : 0.0);

            blGy(i,j) = dy>0.0 ? ((okp ? zp : zc) - (okm ? zm : zc))/dy : 0.0;
        }
    }
}

// grain velocity (ub,vb) and excess Shields number te of a grain d in the column (ic,jc), true if mobile
bool CPM::bedload_grain(lexer *p, int ic, int jc, double d, double &ub, double &vb, double &te)
{
    const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    const double R = (p->S22 - p->W1)/p->W1;

    d = MAX(d,1.0e-9);

    double rg = (p->S22 - p->W1)*gmag*d;
    double tc = p->Q60*pow(d/p->S20,-0.8);

    // bed slope, at most 0.9 mu_s: beyond the angle of repose the packed bed avalanches by itself,
    // the layer only adds the slope effect to the transport by the flow
    double gx = blGx(ic,jc), gy = p->j_dir==1 ? blGy(ic,jc) : 0.0;
    double gm = sqrt(gx*gx + gy*gy);
    
    if(gm>0.9*mu_s)
    {
        gx *= 0.9*mu_s/gm;
        gy *= 0.9*mu_s/gm;
    }
    
    double tx = blTx(ic,jc)/rg - tc/mu_s*gx;
    double ty = p->j_dir==1 ? blTy(ic,jc)/rg - tc/mu_s*gy : 0.0;
    double tm = sqrt(tx*tx + ty*ty);

    te = tm - tc;

    double sp = p->Q59*sqrt(R*gmag*d)*MAX(0.0, sqrt(tm) - 0.7*sqrt(tc));

    ub = tm>1.0e-12 ? sp*tx/tm : 0.0;
    vb = tm>1.0e-12 ? sp*ty/tm : 0.0;

    return te>0.0 && sp>0.0;
}

// settling velocity, Ferguson & Church (2004) for spheres
double CPM::settling_velocity(lexer *p, double d)
{
    const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    const double R = (p->S22 - p->W1)/p->W1;

    return R*gmag*d*d/(18.0*p->W2 + sqrt(0.75*0.4*R*gmag*d*d*d));
}

// occupancy of the cells (nearest cell, as the grid-limited step) for the jumps of the layer parcels:
// pickup to the bed surface, deposition into the bed and the hops bypass the grid-limited step and must
// keep the cells within their capacity
void CPM::bedload_occupancy(lexer *p)
{
    const double vpar = P.ParcelFactor*Vp;

    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    Locc(i,j,k) = 0.0;

    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        int ic,jc;
        bedload_column(p,P.X[n],P.Y[n],ic,jc);
        int kc = MAX(0, MIN(p->knoz-1, p->posc_k(P.Z[n])));
        Locc(ic,jc,kc) += vpar;
    }
}

// a layer parcel leaves its cell at (x,y,z) and is placed at the level zl of the column (ic,jc):
// the cell of zl if it has room, else the next cell above with room; returns the new z
double CPM::bedload_place(lexer *p, fdm *a, double x, double y, double z, int ic, int jc, double zl, double d)
{
    const double vpar = P.ParcelFactor*Vp;
    int is,js;

    bedload_column(p,x,y,is,js);
    int ks = MAX(0, MIN(p->knoz-1, p->posc_k(z)));
    Locc(is,js,ks) -= vpar;

    int kc = MAX(0, MIN(p->knoz-1, p->posc_k(zl)));
    const int kc0 = kc;
    const double zl0 = zl;

    while(kc<p->knoz-1)
    {
        int ii=ic, jj=jc, kk=kc;
        
        // never into a solid body: the parcel stays at the bed level (cell over-filled)
        if(p->solidread>0 && a->solid(ii,jj,kk)<0.0)
        {
            kc = kc0;
            zl = zl0;
            break;
        }
        
        double V = p->DXN[ii+marge]*p->DYN[jj+marge]*p->DZN[kk+marge];
        double t0 = p->Q12==2 ? T0e(ii,jj,kk) : theta_0;
        double cap = MAX((t0 + theta_max - theta_0)*V, t0*V + vpar*(1.0+1.0e-6));

        if(Locc(ii,jj,kk) + vpar <= cap)
        break;

        ++kc;
        zl = MAX(zl, p->ZN[kc+marge] + 0.5*MIN(d, p->DZN[kc+marge]));
    }

    Locc(ic,jc,kc) += vpar;

    return zl;
}

// pickup from the exposed top layer and release into the suspension, once per sub-step
void CPM::bedload_exchange(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, double dt)
{
    const int nj = p->knoy+2;
    const int nc = (p->knox+2)*nj;
    const double vpar = P.ParcelFactor*Vp;
    const bool hasexpo = expo.size()==size_t(P.index);

    std::vector<std::vector<int>> cand(nc), mov(nc);
    std::vector<std::vector<double>> rate(nc);

    int ic,jc;
    double ub,vb,te;

    bedload_occupancy(p);

    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        bedload_column(p,P.X[n],P.Y[n],ic,jc);
        int c = (ic+1)*nj + jc+1;

        if(P.Hop[n]>0.0)
        {
            mov[c].push_back(n);
            continue;
        }

        if(!hasexpo || expo[n]<=0.0)
        continue;

        if(!bedload_grain(p,ic,jc,P.D[n],ub,vb,te))
        continue;

        // share of the bed area times V_m* u_b/L
        double r = expo[n]*vpar/(theta_0*P.D[n]) * (te*P.D[n]/mu_s) * sqrt(ub*ub + vb*vb)/(p->Q62*P.D[n]);

        cand[c].push_back(n);
        rate[c].push_back(r);
    }

    std::uniform_real_distribution<double> uni(0.0,1.0);

    // suspension
    const double ws = settling_velocity(p,p->S20);
    const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    const double R = (p->S22 - p->W1)/p->W1;
    const double Dst = p->S20*pow(R*gmag/(p->W2*p->W2), 1.0/3.0);

    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        int c = (i+1)*nj + j+1;
        auto &cn = cand[c];
        auto &rn = rate[c];

        double rsum = 0.0;
        for(double r : rn)
        rsum += r;

        // morphological factor Q 65: the bed evolves Q 65 times faster than the flow (bedload only)
        blC(i,j) += p->Q65*rsum*dt;

        while(blC(i,j)>=vpar && !cn.empty())
        {
            // weighted choice among the exposed parcels
            double tot = 0.0;
            for(double r : rn)
            tot += r;

            double u = uni(rng)*tot;
            size_t q=0;

            for(q=0; q+1<cn.size(); ++q)
            {
                u -= rn[q];
                if(u<=0.0)
                break;
            }

            n = cn[q];

            // hop length, exponential with the mean Q 62 d
            P.Hop[n] = -p->Q62*P.D[n]*log(MAX(1.0-uni(rng), 1.0e-12));
            P.Z[n] = bedload_place(p,a,P.X[n],P.Y[n],P.Z[n],i,j,zbl(s,i,j)+0.5*P.D[n],P.D[n]);
            P.U[n] = P.V[n] = P.W[n] = 0.0;
            mov[c].push_back(n);
            ++bl_npick;

            cn.erase(cn.begin()+q);
            rn.erase(rn.begin()+q);

            blC(i,j) -= vpar;
        }

        // no grains to pick up: no credit beyond one parcel
        blC(i,j) = MIN(blC(i,j), vpar);

        // below the threshold the remainder decays (1 s)
        if(rsum<=0.0)
        blC(i,j) *= MAX(0.0, 1.0-dt/1.0);

        // release into the suspension
        if(p->Q58==2)
        {
            double tm = sqrt(blTx(i,j)*blTx(i,j) + blTy(i,j)*blTy(i,j));
            double ustar = sqrt(tm/p->W1);
            double T = tm/((p->S22-p->W1)*gmag*p->S20*p->Q60) - 1.0;

            // onset of suspension after van Rijn (1984): u* > 4 w_s/D* (1 < D* <= 10), 0.4 w_s (D* > 10)
            double ucs = (Dst<=10.0 ? 4.0/MAX(Dst,1.0) : 0.4)*ws;
            
            if(ustar>ucs && T>0.0)
            {
                double ha = 0.5*blH(i,j);
                double ca = MIN(0.05, 0.015*p->S20*pow(T,1.5)/(ha*pow(Dst,0.3)));

                // morphological factor Q 65: the suspended parcels act on the bed only (no feedback on the
                // flow with S 10 1), so the release is scaled as the pickup, the deposition follows
                blCs(i,j) += p->Q65*ws*ca*p->DXN[IP]*p->DYN[JP]*dt;

                auto &mn = mov[c];

                while(blCs(i,j)>=vpar && !mn.empty())
                {
                    size_t q = MIN(mn.size()-1, size_t(uni(rng)*mn.size()));
                    n = mn[q];

                    double zr = zbl(s,i,j) + ha;

                    P.Hop[n] = 0.0;
                    P.Z[n] = bedload_place(p,a,P.X[n],P.Y[n],P.Z[n],i,j,zr,P.D[n]);
                    P.U[n] = p->ccipol1c(a->u,P.X[n],P.Y[n],zr);
                    P.V[n] = p->j_dir==1 ? p->ccipol2c(a->v,P.X[n],P.Y[n],zr) : 0.0;
                    P.W[n] = 0.0;
                    ++bl_nsus;

                    mn.erase(mn.begin()+q);
                    blCs(i,j) -= vpar;
                }

                blCs(i,j) = MIN(blCs(i,j), vpar);
            }
        }
    }
}

// a parcel of the bedload layer: hop along the bed with the grain velocity, deposition at the end of the hop
// or below the threshold; tentative position and velocity in XRK1, URK1 (Test = -1 marks the parcel as done);
// a hop into a solid body (solid level set < d/2 at the target) ends where the parcel is
void CPM::bedload_move(lexer *p, fdm *a, sediment_fdm *s, int q, double dt)
{
    int ic,jc;
    double ub,vb,te;

    bedload_column(p,P.X[q],P.Y[q],ic,jc);

    bool mob = bedload_grain(p,ic,jc,P.D[q],ub,vb,te);
    double sp = sqrt(ub*ub + vb*vb);
    double f = 1.0;

    bool dep = !mob || P.Hop[q] <= sp*dt;

    if(mob && dep && sp*dt>1.0e-20)
    f = P.Hop[q]/(sp*dt);

    if(!mob)
    f = 0.0;

    P.XRK1[q] = P.X[q] + f*dt*ub;
    P.YRK1[q] = p->j_dir==1 ? P.Y[q] + f*dt*vb : P.Y[q];

    bedload_column(p,P.XRK1[q],P.YRK1[q],ic,jc);
    double zb = zbl(s,ic,jc);

    if(p->solidread>0 && p->ccipol4_b(a->solid,P.XRK1[q],P.YRK1[q],zb+0.5*P.D[q]) < 0.5*P.D[q])
    {
        P.XRK1[q] = P.X[q];
        P.YRK1[q] = P.Y[q];
        bedload_column(p,P.X[q],P.Y[q],ic,jc);
        zb = zbl(s,ic,jc);
        dep = true;
    }

    P.Test[q] = -1.0;

    if(dep)
    {
        // deposited on the bed, just inside the bed surface (or the first cell above with room)
        double zl = zb - 0.5*P.D[q];

        if(p->nb5<0)
        zl = MAX(zl, p->ZN[0+marge] + 0.5*P.D[q]);

        P.Hop[q] = 0.0;
        P.ZRK1[q] = bedload_place(p,a,P.X[q],P.Y[q],P.Z[q],ic,jc,zl,P.D[q]);
        P.URK1[q] = P.VRK1[q] = P.WRK1[q] = 0.0;
        ++bl_ndep;
    }
    else
    {
        P.Hop[q] -= sp*dt;
        P.ZRK1[q] = bedload_place(p,a,P.X[q],P.Y[q],P.Z[q],ic,jc,zb + 0.5*P.D[q],P.D[q]);
        P.URK1[q] = ub;
        P.VRK1[q] = vb;
        P.WRK1[q] = 0.0;
    }
}

// a grain at rest in the bed or on the bed surface (within one cell above the bed level): with the
// bedload layer it is moved by the layer only (no resolved drag, no turbulent dispersion)
bool CPM::bedload_rest(lexer *p, fdm *a, int q)
{
    double topo = ptopo(p,a,P.X[q],P.Y[q],P.Z[q]);

    if(topo<0.0)
    return true;

    i = p->posc_i(P.X[q]);
    j = p->posc_j(P.Y[q]);
    k = p->posc_k(P.Z[q]);

    double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
    double sp = sqrt(P.U[q]*P.U[q] + P.V[q]*P.V[q] + P.W[q]*P.W[q]);

    return topo<h && sp<0.1*settling_velocity(p,P.D[q]);
}

// bed level of a column for the layer: the actual level from the parcels (the fluid sees it relaxed in time, Q 63)
double CPM::zbl(sediment_fdm *s, int ic, int jc)
{
    return zbl_ok ? blZb(ic,jc) : s->bedzh(ic,jc);
}
