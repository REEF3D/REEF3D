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
// when the largest front-face slope falls below A 360 deg. Its front face and the band of A 502 L_k behind the crest
// are the seed of the viscosity, with the wavenumber k = 2 pi / L of the crest (L trough to trough).
//
// Viscosity (breaking_kennedy with the seed): the eddy viscosity of Kennedy et al. (2000) as in A 550 3, continuous in
// deta/dt, with three changes:
//  - an event starts in a seeded column or, as in A 550 3, where deta/dt > A 554 c. In a seeded column its initial
//    threshold is the larger of A 554 c and deta/dt of the column at the start, so B starts from 0, and the event
//    ages also below the threshold: the threshold relaxes to eta_t^F = A 501 c over T* from the onset of the crest,
//    and the viscosity grows continuously (starting with eta_t^F at once radiated short waves offshore in Ting &
//    Kirby). The event moves with the front through the age carried from the neighbours, as in A 550 3. With the seed
//    as the only onset, the noise of the crest speed C moved the onset from wave to wave by up to 2.6 m in Ting &
//    Kirby and the wave heights oscillated from gauge to gauge; the tracking adds the onsets that A 550 3 misses
//  - the depth h of the formula is the length L_k = tanh(k h)/k of the wave of the event (h in shallow water, 1/k in
//    deep water), c = sqrt(g L_k) the linear phase speed, T* = A 503 sqrt(L_k/g): nu_b = B delta_b^2 L_k deta/dt
//    over a band of A 502 L_k (also the seed behind the crest)
//  - nu_b f(z) with the depth profile of the orbital velocity f(z) = cosh(k (z - z_b))/cosh(k h)
//    (nhflow_breaking_profile.h), in all layers on all faces of the implicit diffusion.
// The viscosity of the tracked crests themselves (energy matched, as FNPF A 350 4) switched on and off with the state
// of every crest point and radiated short waves in Ting & Kirby; the seed keeps the onset of the tracking and the
// smooth viscosity of A 550 3.
// No breaking in the relaxation zones (ioflow::relax_boxes). Crest log with P 59 1.
//
// A 550 5: A 550 4 plus the deep-water breakers of FNPF A 350 5 (breaking_crest3d model 5, the same keys A 367, A 369,
// A 404, A 407-A 409, A 411): crests in deep and intermediate water (C^2/(g h) < 0.6, 2 pi h/L > 1.4, ramps to 0.9 and
// 0.8) break at B > A 359 - A 369/sqrt(N) (the correction fades out between N = 120 and 200 cells per wavelength), with
// the breaker scales at onset (L_b, H_b, c_b, T_b), the strength eps = b c_b^5/g, b = 0.4 (S - 0.08)^2.5 of the
// breaker slope S = k H/2 (Romero et al. 2012; A 367 2: b(Gamma) of Derakhti et al. 2018), b <= A 407, for A 408 T_b
// and at most A 409 of g H_b^2 L_b/8. eps is dissipated by the breaking viscosity alone, nu = eps / (dissipation per
// unit viscosity of the resolved velocities, nu f(z) on the painted interval), with the depth profile of
// nhflow_vbz; a deep-water breaker does not seed the Kennedy viscosity. The bores of shallow water as A 550 4.
// Energy loss of the focused packets of Rapp & Melville (L_c/100, 16 layers, explicit; A 550 4 and the LES of Derakhti &
// Kirby 2016 in brackets): RN 0.2 % (0), RS1 15.7 % (8.0 %, 14.0 %), RP1 27.3 % (26.2 %, 28.9 %), RP2 31.3 % (30.1 %,
// 31.7 %). The semi-implicit step loses about as much without a breaking model (RS1 10.3 %, RP1 26.5 %, RP2 29.5 %):
// its crests do not reach the onset on RS1 (B <= 0.48) and the viscosity changes RP1 and RP2 by 1-2 points.

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
        
        // A 550 5: the deep-water breakers of FNPF A 350 5 (fnpf_breaking_barthelemy.cpp, breaking_crest3d.h) with the
        // same keys: onset A 359 - A 369/sqrt(N), b(S) of Romero et al. (2012) (A 367 1) or b(Gamma) of D18 (A 367 2),
        // b <= A 407, duration A 408 T_b, energy <= A 409 of the breaker, limiting slope A 404, window A 411. The
        // strength eps is dissipated by the breaking viscosity alone (A 412 is not used): NHFLOW has the velocity
        // gradients of the roller, so the viscosity is matched to eps directly
        if(p->A550==5)
        {
            bt->model = 5;
            bt->BonN = MAX(p->A369,0.0);
            bt->deep = (p->A367==1) ? 2 : ((p->A367==2) ? 1 : 0);
            bt->bmax = MAX(p->A407,0.0);
            bt->tau = MAX(p->A408,0.0);
            bt->win = MAX(p->A411,1.0e-6);
            bt->Slim = (p->A404>0.0) ? p->A404 : 1.0e20;
            bt->fvisc = 1.0;
            
            // candidates (rays, log) below the lowest onset threshold (N = 30) and its B window
            bt->Bcand = MIN(bt->Bcand, 0.9*(bt->Bon - bt->BonN/sqrt(30.0) - bt->win));
            
            if(p->mpirank==0)
            cout<<"A 550 5: deep-water breakers of A 350 5 (b <= "<<bt->bmax<<", duration "<<bt->tau<<" T_b, energy <= "
                <<p->A409<<"), dissipated by the breaking viscosity; bores as A 550 4"<<endl;
        }
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

    bt->detect(p,pgc,eta_n,bt_us,bt_vs,d->WL,brkzone_built ? brkzone.data() : nullptr,t);

    // the seed: front faces and bands of the breaking crest points, with their wavenumber
    // painted behind the crest: A 502 L_k, L_k = tanh(k h)/k (h in shallow water; A 502 h reached over the trough
    // behind the crest in deep water)
    auto sback = [&](const bc3_point &R) { return p->A502*((R.kw>0.0) ? tanh(R.kw*R.wl)/R.kw : R.wl); };
    auto integrand = [&](const bc3_point &R, double x, double y) { return 1.0; };
    auto epsf = [&](const bc3_point &R) { return 1.0; };

    SLICELOOP4
    {
        bt_vbs(i,j) = 0.0;
        brkkseed(i,j) = 0.0;
        brkseed(i,j) = 0;
    }

    if(p->A550!=5)
    bt->rays(p,pgc,eta_n,d->WL,t,sback,integrand,epsf,bt_vbs,brkseed.data(),&brkkseed);
    
    // A 550 5: the deep-water breakers (mech 1) get the viscosity nu = eps / (dissipation per unit viscosity) in bt_vbd,
    // with the depth profile of nhflow_vbz (wavenumber k of the crest); the bores (mech 0) seed the Kennedy viscosity as
    // with A 550 4 (their viscosity bt_vbs is not used)
    if(p->A550==5)
    {
        const double g = fabs(p->W22);
        const bool si = (p->A510>=11);
        
        // dissipation of nu f(sig) per unit viscosity, crest length and ray length at (x, y): the column of the nearest
        // cell, int f (|grad u|^2 + |grad v|^2 + |grad w|^2) dz (sigma layers, horizontal derivatives along the layers),
        // the normal stresses with the factor 2 of the explicit implicit diffusion (nhflow_idiff); the semi-implicit
        // horizontal and vertical diffusion is of the Laplacian form (factor 1)
        const double fn = si ? 1.0 : 2.0;
        
        auto integrand5 = [&](const bc3_point &R, double x, double y) -> double
        {
            int ii = bt->cell_i(p,x), jj = bt->cell_j(p,y);
            
            // the nearest cell centre
            if(ii>=0 && ii+1<p->knox && fabs(p->XP[ii+1+marge]-x)<fabs(p->XP[ii+marge]-x))
            ++ii;
            
            if(p->j_dir==1 && jj>=0 && jj+1<p->knoy && fabs(p->YP[jj+1+marge]-y)<fabs(p->YP[jj+marge]-y))
            ++jj;
            
            if(ii<0 || ii>=p->knox || (p->j_dir==1 && (jj<0 || jj>=p->knoy)))
            return 0.0;
            
            const int jc = (p->j_dir==1) ? jj : 0;
            
            if(p->wet[(ii-p->imin)*p->jmax + (jc-p->jmin)]==0)
            return 0.0;
            
            const double h = d->WL(ii,jc);
            
            if(h<=p->A544)
            return 0.0;
            
            const double kh = R.kw*h;
            const double dxx = p->XP[ii+1+marge] - p->XP[ii-1+marge];
            const double dyy = (p->j_dir==1) ? p->YP[jc+1+marge] - p->YP[jc-1+marge] : 1.0;
            double D = 0.0;
            
            for(int kk=0; kk<p->knoz; ++kk)
            {
                const double sg = p->ZP[kk+marge];
                const double f = (kh>0.0) ? (exp(kh*(sg-1.0)) + exp(-kh*(sg+1.0)))/(1.0 + exp(-2.0*kh)) : 1.0;
                const int km = MAX(kk-1,0), kp = MIN(kk+1,p->knoz-1);
                const double dz = h*(p->ZP[kp+marge] - p->ZP[km+marge]);
                
                const double *F3[3] = {d->U, d->V, d->W};
                double phi = 0.0;
                
                for(int c=0; c<3; ++c)
                {
                    if(c==1 && p->j_dir==0)
                    continue;
                    
                    const double *A = F3[c];
                    const double ax = (A[idx(ii+1,jc,kk)] - A[idx(ii-1,jc,kk)])/dxx;
                    const double ay = (p->j_dir==1) ? (A[idx(ii,jc+1,kk)] - A[idx(ii,jc-1,kk)])/dyy : 0.0;
                    const double az = (kp>km && dz>0.0) ? (A[idx(ii,jc,kp)] - A[idx(ii,jc,km)])/dz : 0.0;
                    
                    phi += ((c==0) ? fn : 1.0)*ax*ax + ((c==1) ? fn : 1.0)*ay*ay + ((c==2) ? fn : 1.0)*az*az;
                }
                
                D += f*phi*h*p->DZN[kk+marge];
            }
            
            return D;
        };
        
        // not used (fvisc 1): breaking_crest3d decides the mechanism with a pressure integrand
        const std::function<double(const bc3_point&, double, double)> integrand_pa = [&](const bc3_point &R, double x, double y) { return 0.0; };
        
        // strength per unit crest length and density: a deep-water breaker b c_b^5/g (D18) over A 408 T_b, at most
        // A 409 of its energy g H_b^2 L_b/8; a bore (only its painted interval, the seed, is used) as FNPF
        auto epsf5 = [&](const bc3_point &R) -> double
        {
            if(R.mech==1)
            {
                const double eps5 = R.b*pow(R.cb,5.0)/g;
                const double ecap = p->A409*g*R.Hb*R.Hb*R.Lb/8.0/MAX(p->A408*R.Tb,1.0e-6);
                
                return MIN(eps5, ecap);
            }
            
            return 1.0;
        };
        
        SLICELOOP4
        {
            bt_vbd(i,j) = 0.0;
            bt_kwd(i,j) = 0.0;
            bt_pa(i,j) = 0.0;
        }
        
        bt->rays(p,pgc,eta_n,d->WL,t,sback,integrand5,epsf5,bt_vbs,brkseed.data(),&brkkseed,&integrand_pa,&bt_pa,&bt_vbd,&bt_kwd);
        
        pgc->gcsl_start4(p,bt_vbd,1);
        pgc->gcsl_start4(p,bt_kwd,1);
    }

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
