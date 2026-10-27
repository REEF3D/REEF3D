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

// A 550 4: breaking onset from crest tracking for NHFLOW, with the eddy viscosity of A 550 3 in the deep-water form
//
// Onset: crest tracking in the horizontal plane (breaking_crest3d) at t^n, once per time step. Crest points are the
// maxima of eta along the ridge normal n, the crest speed C comes from the motion of the crest point along n, u_s is
// the horizontal velocity of the top layer extrapolated to the surface. A crest point breaks at B = u_s.n / C > A 359
// (Barthelemy et al. 2018) when it travels at least at half the linear phase speed of the generated waves, and stops
// when the largest front-face slope falls below A 360 deg. Its front face and the band of A 502 h behind the crest
// are the seed of the viscosity, with the wavenumber k = 2 pi / L of the crest (L trough to trough).
//
// Viscosity (breaking_kennedy with the seed): the eddy viscosity of Kennedy et al. (2000) as in A 550 3, continuous in
// deta/dt, with three changes:
//  - an event starts only in a seeded column (no onset from deta/dt > A 554 c). Its initial threshold is the larger of
//    A 554 c and deta/dt of the column at the start, so B starts from 0, and the event ages in the seeded columns also
//    below the threshold: the threshold relaxes to eta_t^F = A 501 c over T* from the onset of the crest, and the
//    viscosity grows continuously (starting with eta_t^F at once radiated short waves offshore in Ting & Kirby).
//    The event moves with the front through the age carried from the neighbours, as in A 550 3
//  - the depth h of the formula is the length L_k = tanh(k h)/k of the wave of the event (h in shallow water, 1/k in
//    deep water), c = sqrt(g L_k) the linear phase speed, T* = A 503 sqrt(L_k/g): nu_b = B delta_b^2 L_k deta/dt
//  - nu_b f(z) with the depth profile of the orbital velocity f(z) = cosh(k (z - z_b))/cosh(k h)
//    (nhflow_breaking_profile.h), in all layers on all faces of the implicit diffusion.
// The viscosity of the tracked crests themselves (energy matched, as FNPF A 350 4) switched on and off with the state
// of every crest point and radiated short waves in Ting & Kirby; the seed keeps the onset of the tracking and the
// smooth viscosity of A 550 3.
// No breaking in the relaxation zones (ioflow::relax_boxes). Crest log with P 59 1.

#include"nhflow_breaking.h"
#include"breaking_crest3d.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>

void nhflow_breaking::breaking_barthelemy(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta_n)
{
    // once per time step, at t^n; the RK stages keep vb
    if(p->count<=bt_count)
    return;

    bt_count = p->count;

    const double t = p->simtime;

    if(bt==nullptr)
    {
        bt = new breaking_crest3d(p);
        bt->Bon = p->A359;
        bt->slope_end = tan(p->A360*PI/180.0);
        bt->deep = p->A367;
        bt->Tdom = (p->B90>0) ? ((p->wTp>0.0) ? p->wTp : p->wT) : 0.0;
    }

    bt->Hmin = 2.0*p->A544;

    // horizontal velocity at the surface: the top layer extrapolated linearly from the two top layers
    const int kt = p->knoz-1;

    auto idx = [&](int ii, int jj, int kk) { return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin; };

    SLICELOOP4
    {
        double u1 = d->U[idx(i,j,kt)];
        double v1 = d->V[idx(i,j,kt)];

        if(kt>0)
        {
            const double z1 = p->ZP[kt+marge], z0 = p->ZP[kt-1+marge];
            const double f = (1.0-z1)/(z1-z0);

            u1 += (u1 - d->U[idx(i,j,kt-1)])*f;
            v1 += (v1 - d->V[idx(i,j,kt-1)])*f;
        }

        bt_us(i,j) = (p->wet[IJ]==1) ? u1 : 0.0;
        bt_vs(i,j) = (p->wet[IJ]==1 && p->j_dir==1) ? v1 : 0.0;
    }

    pgc->gcsl_start4(p,bt_us,1);
    pgc->gcsl_start4(p,bt_vs,1);

    bt->detect(p,pgc,eta_n,bt_us,bt_vs,d->WL,brkzone_built ? brkzone.V : nullptr,t);

    // the seed: front faces and bands of the breaking crest points, with their wavenumber
    auto sback = [&](const bc3_point &R) { return p->A502*R.wl; };
    auto integrand = [&](const bc3_point &R, double x, double y) { return 1.0; };
    auto epsf = [&](const bc3_point &R) { return 1.0; };

    SLICELOOP4
    {
        bt_vbs(i,j) = 0.0;
        brkkseed(i,j) = 0.0;
        brkseed(i,j) = 0;
    }

    bt->rays(p,pgc,eta_n,d->WL,t,sback,integrand,epsf,bt_vbs,brkseed.V,&brkkseed);

    pgc->gcsl_start4int(p,brkseed,50);
    pgc->gcsl_start4(p,brkkseed,1);

    if(p->P59==1)
    {
        if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_Breaking_Log",0777);

        bt->log(p,t,"./REEF3D_NHFLOW_Breaking_Log/REEF3D-NHFLOW-Breaking-Crests.dat");
    }

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking crest points: "<<bt->nbreak<<endl;
}
