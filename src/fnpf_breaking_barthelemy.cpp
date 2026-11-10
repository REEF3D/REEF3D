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
// A352 5: instead of the filter a viscosity floor in the cells of the depth-induced criterion (bt_filter).
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
// call of the step); the RK stages reuse vb. No breaking in the relaxation zones (ioflow::relax_boxes), with A 405 < 1
// only in their outer fraction A 405.
//
// Surf zone settings: A352 0 and a coastline band A341 ~10 (A343 1). The A352 filter (A354 0.6) flattens
// the fronts before B reaches the threshold, and the breaking point is that of A350 1 again; without the
// filter and with the default band A341 2 the swash zone of Ting & Kirby goes unstable after 20-36 s.

// A 350 5: A 350 4 with the breaker scales, strength, duration and energy of steepness-limited breaking from the
// literature (Derakhti, Banner & Kirby 2018, JFM 848, "D18"; Romero et al. 2012; Rapp & Melville 1990; Tian et al.
// 2010) and an energy-matched absorbing surface pressure (Mohanlal et al. 2023), in deep and intermediate water (crests with C^2/(g h) < 0.6 and 2 pi h/L > 1.4,
// ramps to 0.9 and 0.8); bores in shallow water as A 350 4.
//   onset      B > A 359 - A 369/sqrt(N) as A 350 4, but the resolution correction fades out between N = 120 and
//              200 cells per crest wavelength (it was calibrated at N = 30-120; B = 0.85 holds on fine grids)
//   scales     at onset: L_b (trough to trough), H_b, k_b = 2 pi/L_b, c_b = (g tanh(k_b h)/k_b)^1/2 (linear, D18 app. A:
//              not the tracked crest speed, which is 0.7-0.9 c_b at onset), T_b = L_b/c_b
//   strength   eps = b c_b^5/g per unit crest length and density (D18 eq. 1.1), b <= A 407 (0.068, the largest b of
//              the D18 data), with A 367
//              1 (default): b = 0.4 (S - 0.08)^2.5, the laboratory fit of Romero, Melville & Kleiss (2012, JPO 42)
//                (Drazen et al. 2008: b ~ S^2.5), S = k H/2 of the breaker (H crest to trough ahead, k = 2 pi/L), the
//                largest S while the B window of 2 is open (until B > Bth + A 411, at most 0.25 T_b): geometric, so independent of how well the grid resolves
//                the growth of B. Packets of Tian et al. (2010) and Rapp & Melville (1990) at L/40 against the LES
//                losses of Derakhti & Kirby (2016): mean |error| 3.6 points (A 350 4: 3.1), no false onset
//              2: Gamma = T_b dB/dt from a least-squares fit of B(t) over |B - Bth| < A 411 (0.03, D18 sec. 3) of the
//                crest's own samples (every step, before and after onset), x A 359/Bth on coarse grids,
//                b = 0.034 (Gamma - 0.30)^2.5 (D18 eq. 3.1). At L/40 Gamma comes out at 0.2-0.5 for the Rapp &
//                Melville packets and the losses are too low (RS1 2.5 % against 14 %): Gamma needs a resolved B(t)
//   duration   A 408 T_b (0.75, the averaging interval of the b fit, D18 app. A), then the breaker ends; the crest
//              breaks again only after its B has fallen below the window (one event per steepening)
//   limit      A 404 (0.44): the limiting slope k H/2 of a Stokes wave (H/L = 0.142). A crest steeper than 1.25 A 404
//              (0.55, the largest front-face slopes of Toffoli et al. 2010) starts breaking whatever its B, a breaker
//              steeper than A 404 does not end and is absorbed with nua = c_b (on grids that do not resolve the
//              growth of B, B can reach the threshold only after the crest has passed the limiting steepness, and the
//              crest then grows until the run stops; 0: off)
//   energy     at most A 409 (0.25) of g H_b^2 L_b/8 per event (Rapp & Melville 1990: up to 25 % for plungers,
//              Tian et al. 2010: 8-25 %): eps <= A 409 g H_b^2 L_b / (8 A 408 T_b)
//   dissipation A 412: 1 the eddy viscosity on the front face as A 350 4, matched to eps (Tian et al. 2012,
//              Seiffert & Ducrozet 2018); 3 (default) half of eps each; 2 absorbing surface pressure pa = nua S eta_t on the front face (crest to
//              trough ahead, sin^2 ramps) in the dynamic FSBC (Guignard & Grilli 2001, Mohanlal et al. 2023),
//              nua = eps / int S eta_t^2 ds with eta_t of the step, nua <= c_b and <= 0.25 dx/dt (explicit term);
//              The pressure alone does not keep the front face of a breaking potential-flow
//              crest regular: in the 8.6 km JONSWAP run (dx 1.5 m) front-face slopes reached 2 and the run stopped at
//              t = 4431 s; the viscosity damps the sharp front (rate 2 nu k^2 per FSBC). Packets at L/40 (loss against
//              the LES, TP1/TP2/TP3/RS1/RP1/RP2): A 412 1 mean |error| 3.7 points (RP2 -8.1), 2: 3.5, 3: 3.2.
//   log        P 59 1 appends Gamma nwin c_b T_b L_b H_b mech nua dpa (mech 1: deep-water breaker (A 412), 0: bore)

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
    if(bt_zone_built || (p->A350!=4 && p->A350!=5) || pflow==nullptr)
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

    // A 405 < 1: only the outer fraction A 405 of each zone (from its end at the domain boundary) stays free of crest
    // points. The target weight 1 - (exp(x^3.5) - 1)/(e - 1) of the relaxation is below 5 % over the inner half of
    // a zone (x < 0.5), and a steep sea develops there as in the free domain: in a 2D JONSWAP sea (Hs 15.6 m, Tp 15 s,
    // 900 m generation zone) crests grew to 20 m with front slopes of 3 in the last 80 m of the zone, where A 350 4
    // did not act, and blew up at its end (patch 0035)
    if(p->A405<1.0)
    {
        const double fr = MAX(p->A405,0.0);
        const double tol = 1.0e-6*(p->global_xmax-p->global_xmin);
        
        for(size_t q=0; q+3<b.size(); q+=4)
        {
            const bool atxs = b[q]<=p->global_xmin+tol, atxe = b[q+1]>=p->global_xmax-tol;
            const double lx = b[q+1]-b[q];
            
            // a zone along x (at x-min or x-max): shrink it from its inner side
            if(atxs && !atxe)
            b[q+1] = b[q] + fr*lx;
            else
            if(atxe && !atxs)
            b[q] = b[q+1] - fr*lx;
            
            if(p->j_dir==1)
            {
            const double toly = 1.0e-6*MAX(p->global_ymax-p->global_ymin,1.0e-12);
            const bool atys = b[q+2]<=p->global_ymin+toly, atye = b[q+3]>=p->global_ymax-toly;
            const double ly = b[q+3]-b[q+2];
            
            // a zone along y spanning x
            if(atxs && atxe)
            {
            if(atys && !atye)
            b[q+3] = b[q+2] + fr*ly;
            else
            if(atye && !atys)
            b[q+2] = b[q+3] - fr*ly;
            }
            }
        }
        
        if(p->mpirank==0)
        cout<<"A 350 4: crest points also in the inner "<<100.0*(1.0-fr)<<" % of the relaxation zones (A 405)"<<endl;
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

        // A 350 5: Derakhti et al. (2018) strength and duration, absorbing pressure (breaking_crest3d.h)
        if(p->A350==5)
        {
            bt->model = 5;
            // A 367: 1 the slope law of Romero et al. (2012) (default), 2 the Gamma law of D18, 0 off
            bt->deep = (p->A367==1) ? 2 : ((p->A367==2) ? 1 : 0);
            bt->bmax = MAX(p->A407,0.0);
            bt->tau = MAX(p->A408,0.0);
            bt->win = MAX(p->A411,1.0e-6);
            bt->Slim = (p->A404>0.0) ? p->A404 : 1.0e20;
            // A 412: dissipation of the deep-water breakers, 1 viscosity, 2 absorbing pressure, 3 half each (default)
            bt->fvisc = (p->A412==1) ? 1.0 : ((p->A412==2) ? 0.0 : 0.5);
            // the window below the lowest onset threshold is sampled as well
            bt->Bcand = MIN(bt->Bcand, 0.9*(bt->Bon - bt->BonN/sqrt(30.0) - bt->win));
        }
    }

    // A 350 5: eta_t of the step (the absorbing pressure is matched to eps with it), the explicit limit of the pressure
    if(p->A350==5)
    {
        const double rdt = 1.0/MAX(alpha*p->dt,1.0e-20);

        SLICELOOP4
        bt_et(i,j) = (eta(i,j)-eta_n(i,j))*rdt;

        pgc->gcsl_start4(p,bt_et,1);

        bt->nuacap = 0.25*p->DXM/MAX(p->dt,1.0e-20);
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

    // A 350 5: dissipation of the absorbing pressure per unit coefficient, int S eta_t^2 ds
    const std::function<double(const bc3_point&, double, double)> integrand_pa = [&](const bc3_point &R, double x, double y)
    {
        const double e = bt->ip(p,bt_et,x,y);

        return std::isnan(e) ? 0.0 : e*e;
    };

    auto epsf = [&](const bc3_point &R)
    {
        // A 350 5, a deep-water breaker (absorbing pressure): eps = b rho c_b^5 / g over tau T_b (D18), the energy of
        // the event at most A 409 of the breaker's g H_b^2 L_b / 8
        if(R.mech==1)
        {
            const double eps5 = R.b*pow(R.cb,5.0)/g;
            const double ecap = p->A409*g*R.Hb*R.Hb*R.Lb/8.0/MAX(p->A408*R.Tb,1.0e-6);

            return MIN(eps5, ecap);
        }

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
        bt_pa(i,j) = 0.0;
    }

    if(p->A350==5)
    bt->rays(p,pgc,c->eta,c->WL,t,sback,integrand,epsf,c->vb,c->breaking.V,nullptr,&integrand_pa,&bt_pa);
    else
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
// but only shoreward of the deepest breaking onset so far (still-water depth < bt_db, not with A352 5): in the surf and swash
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

            if(p->wet[IJ]==1 && (p->wd-c->bed(i,j)<bt_db || p->A352==5) && (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A354*sqrt(9.81*c->WL(i,j)))
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

    // A352 5: no filter; the breaking viscosity is raised to at least nu = dx c/(4 C), C = 0.5, in the cells of the
    // depth-induced criterion (flag 2) and +-2 cells around them, c the phase speed of the 2 dx wave in the local
    // depth: a 2 dx spike decays by e while it travels half a cell, longer waves much less (about with the square
    // of their length). The maximum with the breaking and coastline viscosity, applied by fnpf_fsfbc::damping.
    // The depth-induced cells are not limited to the still-water depth of the deepest onset (bt_db): at a steep
    // shore the crest tracking sees the onset late and shallow, and the swash further out went unstable
    // (validation 09). C 0.5 and no depth limit: 29 of 30 runs stable on two 2D steep-shore cases (A352 3: 27).
    if(p->A352==5)
    {
        const double g = fabs(p->W22);
        const double C = 0.5;

        SLICELOOP4
        bx(i,j) = (c->breaking(i,j)>=2 && p->wet[IJ]==1) ? 1 : 0;

        pgc->gcsl_start4int(p,bx,50);

        SLICELOOP4
        if(p->wet[IJ]==1)
        {
            int hit=0;

            for(int qi=-2; qi<=2; ++qi)
            for(int qj=-2*p->j_dir; qj<=2*p->j_dir; ++qj)
            if(bx(i+qi,j+qj)==1)
            hit=1;

            if(hit==1)
            {
                const double wl = MAX(c->WL(i,j), c->wd_criterion);
                const double dx = (p->j_dir==1) ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];
                const double k = PI/dx;
                const double cp = sqrt(g*tanh(k*wl)/k);

                c->vb(i,j) = MAX(c->vb(i,j), dx*cp/(4.0*C));
            }
        }

        pgc->gcsl_start4(p,c->vb,1);
        return;
    }

    SLICELOOP4
    if(c->breaking(i,j)>=1 && p->wet[IJ]==1)
    {
        filter(p,c,pgc,eta,p->A361,p->A362);
        filter(p,c,pgc,Fifsf,p->A361,p->A362);
    }
}

// A 350 5: absorbing surface pressure of the deep-water breakers in the dynamic FSBC (Guignard & Grilli 2001,
// Mohanlal et al. 2023), Fi_t -= pa with pa = bt_pa eta_t: its work on the surface is -int bt_pa eta_t^2 dx <= 0
// (eta_t = phi_n sqrt(1 + |grad eta|^2)); bt_pa is set once per step from the state at t^n, ek: eta_t of the stage
void fnpf_breaking::absorb(lexer *p, fdm_fnpf *c, slice &ek, slice &fk)
{
    if(p->A350!=5)
    return;

    SLICELOOP4
    fk(i,j) -= bt_pa(i,j)*ek(i,j);
}
