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

#include"seastate_f.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"slice4.h"
#include"seastate_obstacle.h"
#include"seastate_amr.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<iostream>

/*--------------------------------------------------------------------
Phase-decoupled diffraction (A 718), Holthuijsen, Herman and Booij
(2003), as SWAN DIFFRACTION (structured grid, SWAN DIFPAR):

  delta = div(cg/k grad(sqrt E)) / (k cg sqrt E),   Ca = sqrt(1 + delta)  (1 where delta <= -1)

A 718 1 (SWAN): E the total energy of the cell, k and cg its energy-
weighted means, one Ca for all frequencies. A 718 2: per frequency with
the energy, k and cg of the frequency. E is smoothed n times with the
active neighbours, E - 0.2 sum (E - E_neighbour) (SWAN smpar 0.2, smnum
n) before the square root; n = A 719, or (A 719 0) n = 0.4 (L/dx)^2 from
the mean wavelength L = 2 pi/k of the domain and the smallest cell size
(smoothing over a fixed fraction of the wavelength: on finer grids the
curvature of sqrt E is noisier and Ca starts to oscillate), 1 to 400. Ca = 1 where E is below 10^-6 of
its maximum. With mesh refinement (G 1, seastate_amr::diffraction) every patch smooths its own energy
with 0.4 (L/dx)^2 steps of its cell size (A 719 n: n 4^level), the ring around its interior held at the
smoothed energy and Ca of the next coarser grid. The solver scales the geographic velocities with Ca and
adds the turning

  c_theta += cg dCa/dn = cg (-sin(theta) dCa/dx + cos(theta) dCa/dy)

to the depth refraction (scaled with Ca), i.e. the rays turn towards
larger Ca (larger K = k Ca), which is the shadow side of a shadow line:
energy spreads behind obstacles and islands.
Gradients central, one-sided next to the edge of the active cells and at
obstacles (A 722), which also bound the smoothing and the divergence (as SWAN).
Computed from the latest spectra before every iteration (stationary)
or step (nonstationary).
--------------------------------------------------------------------*/

void seastate_f::diffraction(lexer *p, ghostcell *pgc)
{
    if(p->A718!=1 && p->A718!=2)
    return;

    const seastate_grid &g = *e->grid;
    slice4 &S = *dfS, &T = *dfT, &KM = *dfK, &CM = *dfC;
    int i,j;

    // number of smoothing steps: A 719, or from the mean wavelength and the smallest cell (the first call
    // with energy in the domain; before, one step)
    int nsmooth = std::max(dfsmooth,1);

    if(dfsmooth<0)
    {
        if(p->A719>0)
        dfsmooth = p->A719;
        else
        {
        double ek = 0.0, et = 0.0, dmin = 1.0e30;
            SLICELOOP4
            if(e->wet(i,j)==1 && e->N->spec(i,j)!=nullptr)
            {
            const float *N = e->N->spec(i,j), *kk = e->kw->spec(i,j);
                for(int l=0; l<g.nsig; ++l)
                {
                double s = 0.0;
                for(int m=0; m<g.ndir; ++m)
                s += double(N[g.bin(l,m)])*g.wth[m];
                s *= g.sig[l]*g.dsig[l];
                et += s;
                ek += s*double(kk[l]);
                }
            dmin = std::min(dmin,std::min(p->DXN[IP],p->DYN[JP]));
            }
        et = pgc->globalsum(et);
        ek = pgc->globalsum(ek);
        dmin = pgc->globalmin(dmin);
        const double L = (ek>0.0 && et>0.0) ? 2.0*3.14159265358979323846*et/ek : 0.0;
        if(L>0.0)
        {
        dfsmooth = std::min(std::max(int(std::ceil(0.4*(L/dmin)*(L/dmin))),1),400);
        dfL = L;
        dfdmin = dmin;
        }
        }

        if(dfsmooth>0 && p->mpirank==0)
        cout<<"SEASTATE diffraction: "<<dfsmooth<<" smoothing steps per iteration"<<(p->A719>0 ? " (A 719)" : " (from the mean wavelength and the cell size, A 719 0)")<<endl;

        nsmooth = std::max(dfsmooth,1);
    }

    auto active = [&](int ii, int jj) {return e->wet(ii,jj)==1 && e->N->spec(ii,jj)!=nullptr;};

    // the neighbour across side s (0 x-, 1 x+, 2 y-, 3 y+) is active and the face not blocked by an obstacle
    auto nb = [&](int ii, int jj, int sd)
    {
        const int ni = ii + (sd==1) - (sd==0), nj = jj + (sd==3) - (sd==2);
        if(!active(ni,nj))
        return false;
        if(pobs==nullptr)
        return true;
        const seastate_obstacle::face *f = (sd==0) ? pobs->east(ii-1,jj) : (sd==1) ? pobs->east(ii,jj) : (sd==2) ? pobs->north(ii,jj-1) : pobs->north(ii,jj);
        return f==nullptr;
    };

    // Ca into T from the energy in S and the wave number and group velocity in KM, CM; the smoothed
    // energy into dfE (mesh refinement: the values next to the patches)
    double smax = 0.0;
    auto ca_field = [&]()
    {
        pgc->gcsl_start4(p,S,50);
        pgc->gcsl_start4(p,KM,50);
        pgc->gcsl_start4(p,CM,50);

        for(int it=0; it<nsmooth; ++it)
        {
            SLICELOOP4
            {
            double t = S(i,j);
                if(active(i,j))
                {
                if(nb(i,j,0)) t -= 0.2*(S(i,j)-S(i-1,j));
                if(nb(i,j,1)) t -= 0.2*(S(i,j)-S(i+1,j));
                if(nb(i,j,2)) t -= 0.2*(S(i,j)-S(i,j-1));
                if(nb(i,j,3)) t -= 0.2*(S(i,j)-S(i,j+1));
                }
            T(i,j) = t;
            }

            SLICELOOP4
            S(i,j) = T(i,j);

            pgc->gcsl_start4(p,S,50);
        }

        if(dfE!=nullptr)
        {
        IMALOOP
        JMALOOP
        (*dfE)(i,j) = S(i,j);
        }

        smax = 0.0;
        SLICELOOP4
        smax = std::max(smax,S(i,j));
        smax = pgc->globalmax(smax);

        SLICELOOP4
        S(i,j) = std::sqrt(std::max(S(i,j),0.0));
        pgc->gcsl_start4(p,S,50);

        auto F = [&](int ii, int jj) {return KM(ii,jj)>0.0 ? CM(ii,jj)/KM(ii,jj) : 0.0;};

        SLICELOOP4
        {
        double ca = 1.0;

            if(active(i,j) && S(i,j)*S(i,j)>1.0e-6*smax && smax>0.0)
            {
            const double k = KM(i,j), cgm = CM(i,j), Fp = F(i,j);
            const double rdx2 = 1.0/(p->DXN[IP]*p->DXN[IP]), rdy2 = 1.0/(p->DYN[JP]*p->DYN[JP]);
            double div = 0.0;

            if(nb(i,j,0)) div += 0.5*(Fp+F(i-1,j))*(S(i-1,j)-S(i,j))*rdx2;
            if(nb(i,j,1)) div += 0.5*(Fp+F(i+1,j))*(S(i+1,j)-S(i,j))*rdx2;
            if(nb(i,j,2)) div += 0.5*(Fp+F(i,j-1))*(S(i,j-1)-S(i,j))*rdy2;
            if(nb(i,j,3)) div += 0.5*(Fp+F(i,j+1))*(S(i,j+1)-S(i,j))*rdy2;

                if(k>0.0 && cgm>0.0)
                {
                const double delta = div/(k*cgm*S(i,j));
                if(delta>-1.0)
                ca = std::sqrt(1.0+delta);
                }
            }

        T(i,j) = ca;
        }

        pgc->gcsl_start4(p,T,50);
    };

    // Ca and its gradient into the stores (frequencies l0..l1)
    auto store = [&](int l0, int l1)
    {
        IMALOOP
        JMALOOP
        {
        float *c = dca->spec(i,j);
            if(c!=nullptr)
            for(int l=l0; l<=l1; ++l)
            c[l] = active(i,j) ? float(T(i,j)) : 1.0f;
        }

        SLICELOOP4
        {
        float *cx = dcax->spec(i,j), *cy = dcay->spec(i,j);

            if(cx==nullptr || cy==nullptr)
            continue;

        double gx = 0.0, gy = 0.0;

            if(active(i,j))
            {
            const bool w = nb(i,j,0), ea = nb(i,j,1), s = nb(i,j,2), n = nb(i,j,3);

            if(w && ea)
            gx = (T(i+1,j)-T(i-1,j))/(p->XP[IP1]-p->XP[IM1]);
            else if(ea)
            gx = (T(i+1,j)-T(i,j))/(p->XP[IP1]-p->XP[IP]);
            else if(w)
            gx = (T(i,j)-T(i-1,j))/(p->XP[IP]-p->XP[IM1]);

            if(s && n)
            gy = (T(i,j+1)-T(i,j-1))/(p->YP[JP1]-p->YP[JM1]);
            else if(n)
            gy = (T(i,j+1)-T(i,j))/(p->YP[JP1]-p->YP[JP]);
            else if(s)
            gy = (T(i,j)-T(i,j-1))/(p->YP[JP]-p->YP[JM1]);
            }

            for(int l=l0; l<=l1; ++l)
            {
            cx[l] = float(gx);
            cy[l] = float(gy);
            }
        }
    };

    if(p->A718==1)
    {
        // total energy, energy-weighted k and cg
        SLICELOOP4
        {
        double et = 0.0, ek = 0.0, ec = 0.0;
            if(active(i,j))
            {
            const float *N = e->N->spec(i,j), *kk = e->kw->spec(i,j), *cc = e->cg->spec(i,j);
                for(int l=0; l<g.nsig; ++l)
                {
                double s = 0.0;
                const float *Nl = N + g.bin(l,0);
                for(int m=0; m<g.ndir; ++m)
                s += double(Nl[m])*g.wth[m];
                const double el = s*g.sig[l]*g.dsig[l]*g.dtheta;
                et += el;
                ek += el*double(kk[l]);
                ec += el*double(cc[l]);
                }
            }
        S(i,j) = et;
        KM(i,j) = et>0.0 ? ek/et : 0.0;
        CM(i,j) = et>0.0 ? ec/et : 0.0;
        }

        ca_field();
        store(0,g.nsig-1);

        if(pamr!=nullptr && pamr->active())
        pamr->diffraction(p,1,0,g.nsig-1,smax,dfL,dfdmin,*dfE,T);
        return;
    }

    for(int l=0; l<g.nsig; ++l)
    {
        SLICELOOP4
        {
        double s = 0.0;
            if(active(i,j))
            {
            const float *N = e->N->spec(i,j) + g.bin(l,0);
            for(int m=0; m<g.ndir; ++m)
            s += double(N[m])*g.wth[m];
            }
        S(i,j) = s;
        KM(i,j) = active(i,j) ? double(e->kw->spec(i,j)[l]) : 0.0;
        CM(i,j) = active(i,j) ? double(e->cg->spec(i,j)[l]) : 0.0;
        }

        ca_field();
        store(l,l);

        if(pamr!=nullptr && pamr->active())
        pamr->diffraction(p,2,l,l,smax,dfL,dfdmin,*dfE,T);
    }
}

// obstacles, structures and coasts (A 722 - A 726): transmission after Goda, d'Angremond and of the porous
// structures from the present spectra, coasts after wetting and drying, before every iteration; also on the
// patches of the mesh refinement
void seastate_f::obstacles(lexer *p, ghostcell*)
{
    if(pobs!=nullptr && pobs->active())
    pobs->update(p,e);

    if(pobs!=nullptr && pamr!=nullptr && pamr->active())
    pamr->obstacles();
}
