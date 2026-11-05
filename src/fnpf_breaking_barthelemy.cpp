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

// A 350 4: crest-tracking breaking model for FNPF
//
// Onset (Barthelemy et al. 2018, JFM 841; Derakhti et al. 2020, shallow water):
//   B = u_s / C > A 359 (0.85), u_s = horizontal surface velocity at the crest along the direction of travel n
//   (= grad Fi_s there, the crest has grad eta . n = 0), C = crest speed along n from tracking the crest over time.
//   Deep and intermediate water (C^2/(g h) < 0.6 and 2 pi h/L > 1.4, ramps to 0.9 and 0.8, h the still-water depth):
//   B > A 359 - A 369/sqrt(N), N = L/dx cells per wavelength of the crest. The threshold 0.85 needs the crest slow-down before breaking resolved (dx ~ L/200); on
//   coarser grids C stays high and u_s low. A 369 1.64 from focused packets without onset at N = 30-120 (Tian et al.
//   2010, Rapp & Melville 1990): midway between the largest B of the non-breaking (TN, RI) and the smallest of the
//   breaking packets (TP1, RS1); 0.59 at N = 40, 0.68 at N = 100, N not below 30. A 369 0: A 359 everywhere.
// Termination: largest front-face slope below tan(A 360 deg) (8 deg, Simon et al. 2018).
// Strength, per unit crest length and density (the larger of the two):
//   bore / hydraulic jump   eps = A 364 * g C H^3 / (4 h_c)   (discharge C h_t times the head loss H^3/(4 h_c h_t))
//   deep water (A 367 1)    eps = b C^5 / g, b = 0.034 (Gamma - 0.30)^2.5 (Derakhti, Banner & Kirby 2018),
//                           Gamma = T_b dB/dt at onset (at most 1.62, b <= 0.068), T_b = L/C
//   H crest to trough ahead, h_c water depth at the crest, L trough to trough.
// Filter: with A352 > 0 the filter of A350 1 (A361, A362) acts on eta and Fifsf in the breaking zones and,
// with A351 1/3, in the cells of the depth-induced criterion of A350 1 (A354) shoreward of the deepest onset.
// Coastline: the coastline / wet-dry front viscosity A346 of A350 1 is applied as well.
// Dissipation: the implicit damping of fnpf_fsfbc::damping with one viscosity per breaking crest point on
// its front face (crest to trough ahead along n, sin^2 ramps at both ends). For a wave of permanent form
// (eta_t = -C eta_n, Fi_t = -C Fi_n) the damping terms dissipate
//   D = nu (fe + ff) C int S eta_n Fi_nn ds,   fe = A 366 (A 357 1, 3), ff = 1 (A 357 1, 2),
// per unit crest length, and nu = eps / D-per-nu matches the dissipation of the zone to eps.
//
// The crests are tracked in the horizontal plane without a preferred direction (breaking_crest3d): crest points
// are the maxima of eta along the ridge normal, every crest point continues the nearest one of the last step, and
// the front faces are rays along n through the subdomains. Tracking and viscosity once per time step at t^n (first
// call of the step); the RK stages reuse vb. No breaking in the relaxation zones (ioflow::relax_boxes).
//
// Surf zone settings: A352 0 and a coastline band A341 ~10 (A343 1). The A352 filter (A354 0.6) flattens
// the fronts before B reaches the threshold, and the breaking point is that of A350 1 again; without the
// filter and with the default band A341 2 the swash zone of Ting & Kirby goes unstable after 20-36 s.

#include"fnpf_breaking.h"
#include"breaking_crest3d.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>

// the relaxation zones of the wave generation and the numerical beach: no crest points there. Without relaxation
// zones (Dirichlet generation B 98 3/4) the B 96 lengths at x-min and, with a beach B 99, at x-max, as before.
void fnpf_breaking::breaking_zones(lexer *p, ioflow *pflow)
{
    if(bt_zone_built || p->A350!=4 || pflow==nullptr)
    return;

    std::vector<double> b;
    pflow->relax_boxes(p,b);

    if(b.empty() && p->B90>0)
    {
        const double big = 1.0e20;
        b.insert(b.end(),{-big, p->global_xmin+p->B96_1, -big, big});

        if(p->B99>0)
        b.insert(b.end(),{p->global_xmax-p->B96_2, big, -big, big});
    }

    SLICELOOP4
    {
    bt_zone(i,j) = 0;

        for(size_t q=0; q+3<b.size(); q+=4)
        if(p->XP[IP]>=b[q] && p->XP[IP]<=b[q+1] && (p->j_dir==0 || (p->YP[JP]>=b[q+2] && p->YP[JP]<=b[q+3])))
        bt_zone(i,j) = 1;
    }

    bt_zone_built = true;
}

void fnpf_breaking::breaking_barthelemy(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    // A352 > 0: the breaking filter of A350 1 on the stage values in the breaking zones (every stage)
    if(p->count<=bt_count)
    {
        bt_filter(p,c,pgc,eta,eta_n,Fifsf,alpha);
        return;
    }

    bt_count = p->count;

    const double g = fabs(p->W22);
    const double t = p->simtime;

    if(bt==nullptr)
    {
        bt = new breaking_crest3d(p);
        bt->Bon = p->A359;
        bt->slope_end = tan(p->A360*PI/180.0);
        bt->deep = p->A367;
        bt->Tdom = (p->B90>0) ? ((p->wTp>0.0) ? p->wTp : p->wT) : 0.0;
        bt->BonN = MAX(p->A369,0.0);

        // candidates (rays, log) below the lowest onset threshold, that of N = 30 cells per wavelength
        if(bt->BonN>0.0)
        bt->Bcand = MIN(bt->Bcand, 0.9*(bt->Bon - bt->BonN/sqrt(30.0)));
    }

    bt->Hmin = 2.0*c->wd_criterion;

    // surface velocity grad Fi_s at t^n
    SLICELOOP4
    {
        bt_us(i,j) = (c->Fifsf(i+1,j)-c->Fifsf(i-1,j))/(p->XP[IP1]-p->XP[IM1]);
        bt_vs(i,j) = (p->j_dir==1) ? (c->Fifsf(i,j+1)-c->Fifsf(i,j-1))/(p->YP[JP1]-p->YP[JM1]) : 0.0;
    }

    pgc->gcsl_start4(p,bt_us,1);
    pgc->gcsl_start4(p,bt_vs,1);

    bt->detect(p,pgc,c->eta,bt_us,bt_vs,c->WL,bt_zone_built ? bt_zone.V : nullptr,t,&c->Fifsf);

    // dissipation of the damping per unit viscosity, crest length and ray length
    const double fe = (p->A357==1 || p->A357==3) ? p->A366 : 0.0;
    const double ff = (p->A357==1 || p->A357==2) ? 1.0 : 0.0;

    auto integrand = [&](const bc3_point &R, double x, double y)
    {
        const double dd = R.d;
        const double em = bt->ip(p,c->eta,x-dd*R.nx,y-dd*R.ny);
        const double ep = bt->ip(p,c->eta,x+dd*R.nx,y+dd*R.ny);
        const double fm = bt->ip(p,c->Fifsf,x-dd*R.nx,y-dd*R.ny);
        const double f0 = bt->ip(p,c->Fifsf,x,y);
        const double fp = bt->ip(p,c->Fifsf,x+dd*R.nx,y+dd*R.ny);

        if(std::isnan(em) || std::isnan(ep) || std::isnan(fm) || std::isnan(f0) || std::isnan(fp))
        return 0.0;

        return (fe+ff)*R.C*((ep-em)/(2.0*dd))*((fp-2.0*f0+fm)/(dd*dd));
    };

    auto epsf = [&](const bc3_point &R)
    {
        const double hc = MAX(R.wl, 1.0e-3);
        // moving bore: discharge through the jump C h_t, head loss H^3/(4 h_c h_t)
        const double eps_bore = p->A364*g*fabs(R.C)*pow(R.H,3.0)/(4.0*hc);
        const double eps_deep = R.b*pow(fabs(R.C),5.0)/g;

        return MAX(eps_bore, eps_deep);
    };

    auto sback = [&](const bc3_point &R) { return 0.0; };

    SLICELOOP4
    {
        c->vb(i,j) = 0.0;
        c->breaking(i,j) = 0;
    }

    bt->rays(p,pgc,c->eta,c->WL,t,sback,integrand,epsf,c->vb,c->breaking.V,nullptr);

    // still-water depth of the deepest onset so far: the surf zone for the filter (A352)
    bt_db = MAX(bt_db, bt->deepest_onset);

    // coastline / wet-dry front viscosity A346 as with A350 1 (A343 1 with A337, A343 2/3 with A332)
    coast_viscosity(p,c,pgc);

    pgc->gcsl_start4(p, c->vb, 1);
    pgc->gcsl_start4int(p, c->breaking, 50);

    SLICELOOP4
    c->breaklog(i,j) = c->breaking(i,j);

    // crest log (P 59 1): crest points with B > 0.5 or breaking
    if(p->P59==1)
    {
        if(p->mpirank==0)
        mkdir("./REEF3D_FNPF_Breaking_Log",0777);

        bt->log(p,t,"./REEF3D_FNPF_Breaking_Log/REEF3D-FNPF-Breaking-Crests.dat");
    }

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking crest points: "<<bt->nbreak<<endl;

    bt_filter(p,c,pgc,eta,eta_n,Fifsf,alpha);

    bt_time = t;
}

// A352 > 0: filter eta and Fifsf in the cells of the breaking zones (as A350 1 with A352 3).
// With A351 1/3 also the cells of the depth-induced criterion of A350 1 (eta_t > A354 sqrt(g h), +-2 cells),
// but only shoreward of the deepest breaking onset so far (still-water depth < bt_db): in the surf and swash
// zone the filter keeps the FNPF surface stable (Ting & Kirby stops at t = 20-36 s without it), offshore of
// the break point it would flatten the crests before B reaches the onset threshold.
void fnpf_breaking::bt_filter(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    if(p->A352<=0)
    return;

    if((p->A351==1 || p->A351==3) && p->count>1)
    {
        SLICELOOP4
        {
            bx(i,j)=0;

            if(p->wet[IJ]==1 && p->wd-c->bed(i,j)<bt_db && (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A354*sqrt(9.81*c->WL(i,j)))
            bx(i,j)=2;
        }

        pgc->gcsl_start4int(p,bx,50);

        SLICELOOP4
        {
            int hit=0;

            for(int q=-2; q<=2; ++q)
            {
                if(bx(i+q,j)==2)
                hit=1;

                if(p->j_dir==1 && bx(i,j+q)==2)
                hit=1;
            }

            if(hit==1)
            c->breaking(i,j) = MAX(c->breaking(i,j), 2);
        }
    }

    SLICELOOP4
    if(c->breaking(i,j)>=1 && p->wet[IJ]==1)
    {
        filter(p,c,pgc,eta,p->A361,p->A362);
        filter(p,c,pgc,Fifsf,p->A361,p->A362);
    }
}
