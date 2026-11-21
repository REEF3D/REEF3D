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

#include"nhflow_breaking.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<vector>
#include"breaking_crest3d.h"

nhflow_breaking::nhflow_breaking(lexer* p, fdm_nhf *d, ghostcell *pgc) : bx(p), by(p), bd(p), brkflag(p), brkzone(p),
                                                                            brkage(p), brkage_n(p), brkvb(p), bt_us(p), bt_vs(p), bt_vbs(p), brkseed(p),
                                                                            brkkseed(p), brkk(p), brkk_n(p), brkr(p), brkr_n(p), brkh(p), brkh_n(p), plk_n(p), plh_n(p), etf(p), etf2(p), bt_vbd(p), bt_kwd(p), bt_pa(p), count_n(-1)
{
    SLICELOOP4
    {
    d->breaking(i,j)=0;
    brkflag(i,j)=0;
    brkzone(i,j)=0;
    brkage(i,j)=0.0;
    brkage_n(i,j)=0.0;
    brkseed(i,j)=0;
    brkkseed(i,j)=0.0;
    brkk(i,j)=0.0;
    brkk_n(i,j)=0.0;
    brkr(i,j)=0.0;
    brkr_n(i,j)=0.0;
    brkh(i,j)=0.0;
    brkh_n(i,j)=0.0;
    plk_n(i,j)=0.0;
    plh_n(i,j)=0.0;
    d->plk(i,j)=0.0;
    d->plh(i,j)=0.0;
    bt_vbd(i,j)=0.0;
    bt_kwd(i,j)=0.0;
    bt_pa(i,j)=0.0;
    }
}

nhflow_breaking::~nhflow_breaking()
{
    delete bt;
}

void nhflow_breaking::breaking(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    if(p->A550==1)
    breaking_baquet(p, d, pgc, eta, eta_n, WL, alpha);
    
    if(p->A550==3)
    breaking_kennedy(p, d, pgc, WL);
    
    if(p->A550==4)
    {
    breaking_barthelemy(p, d, pgc, eta_n);
    breaking_kennedy(p, d, pgc, WL);
    }
    
    // A 550 5: the bores as A 550 4, the deep-water breakers with the energy-matched viscosity of breaking_crest3d
    if(p->A550==5)
    {
    breaking_barthelemy(p, d, pgc, eta_n);
    breaking_kennedy(p, d, pgc, WL);
    breaking_deep(p, d, pgc);
    }
}

void nhflow_breaking::breaking_deep(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    int count=0;
    
    SLICELOOP4
    if(p->wet[IJ]==1 && bt_vbd(i,j)>d->vb(i,j))
    {
    d->vb(i,j) = bt_vbd(i,j);
    d->vbk(i,j) = bt_kwd(i,j);
    d->breaking(i,j) = 1;
    ++count;
    }
    
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4(p,d->vbk,1);
    pgc->gcsl_start4int(p,d->breaking,50);
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0) && count>0)
    cout<<"breaking: deep-water breakers in "<<count<<" columns"<<endl;
}

// A 550 3: eddy viscosity of Kennedy, Chen, Kirby & Dalrymple (2000), J. Waterway Port Coastal Ocean Eng. 126(1), 39-47,
// as in FUNWAVE, with deta/dt of the stage (continuity) and the total depth h of the stage water level:
//
//   nu_b = B delta_b^2 h deta/dt,   delta_b = A 504 (1.2)
//   B = 1 for deta/dt >= 2 eta_t*,  deta/dt / eta_t* - 1 for eta_t* < deta/dt < 2 eta_t*,  0 below
//   eta_t* = eta_t^I + min(t_b/T*, 1) (eta_t^F - eta_t^I),  eta_t^I = A 554 sqrt(g h), eta_t^F = A 501 sqrt(g h),
//   T* = A 503 sqrt(h/g)
//
// t_b is the age of the breaking event in the column: a column with deta/dt > eta_t* takes the largest age of
// itself and its neighbours at the end of the last time step plus dt, so the event moves with the front and its
// threshold relaxes from eta_t^I to eta_t^F while it travels; otherwise the age is 0. The viscosity varies
// continuously with deta/dt (B), so there is no switch of the flow equations and no on/off of a constant viscosity.
//
// nu_b is spread (maximum of the neighbours) over a band of A 502 h around the breaking columns, so that it covers
// the crest behind the front: on the front face alone it smooths the front, which takes away the dissipation of the
// shock-capturing scheme, while the non-hydrostatic pressure keeps the crest, and the crest peaks (Ting & Kirby).
// It acts in all layers on all faces of the implicit diffusion (A 512 2, set by A 550 3); with the horizontal faces
// only, the inner surf zone of Ting & Kirby depended on the number of layers (4: H/h 0.48, 10: 0.59 and noisy).
// No breaking in the relaxation zones (breaking_zones).
void nhflow_breaking::breaking_kennedy(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    // multi-block grid (nhflow_amr::breaking_level0): only the arguments, the parts run on all blocks with the
    // exchanges between the blocks in between
    if(brk_defer)
    {
    brk_WL = &WL;
    return;
    }
    
    kennedy_part(p,d,pgc,WL,0,0);
    kennedy_part(p,d,pgc,WL,1,0);
    
    const int nband = pgc->globalimax(kn_band);
    
    for(int qn=0; qn<nband; ++qn)
    {
    kennedy_part(p,d,pgc,WL,2,qn);
    kennedy_part(p,d,pgc,WL,3,qn);
    }
    
    kennedy_part(p,d,pgc,WL,4,0);
    
    const int count=pgc->globalisum(kn_count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
    
    const int non = pgc->globalisum(on_count);
    const double xon = pgc->globalsum(on_xsum);
    
    if(p->mpirank==0 && non>0 && (p->P12<=1 || p->count%p->P12==0 || p->A542>0.0))
    cout<<"onset: "<<non<<" columns at x = "<<xon/double(non)<<" t = "<<p->simtime<<endl;
    
    if(p->S137>0.0)
    {
    const int npl = pgc->globalisum(pl_count);
    const double xpl = pgc->globalsum(pl_xsum);
    
    if(p->mpirank==0 && npl>0)
    cout<<"plunge: "<<npl<<" columns at x = "<<xpl/double(npl)<<" t = "<<p->simtime<<endl;
    }
}

// the parts of breaking_kennedy: 0 the ages of the last step, 1 the viscosity of the columns (kn_count) and the width
// of the band (kn_band cells), 2 the ghost cells before the pass qn of the band, 3 the pass, 4 the flags
void nhflow_breaking::kennedy_part(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, int part, int qn)
{
    const double g = fabs(p->W22);
    
    if(part==0)
    {
    // new time step: the ages (and with A 550 4 the wavenumbers) at the end of the last step
    if(p->count>count_n)
    {
    SLICELOOP4
    {
    brkage_n(i,j) = brkage(i,j);
    brkk_n(i,j) = brkk(i,j);
    brkr_n(i,j) = brkr(i,j);
    brkh_n(i,j) = brkh(i,j);
    
    // S 142: the impact turbulence decays over the time S 142 sqrt(H_b/g) (0: only in the impact window)
    const double td = p->S142*sqrt(MAX(d->plh(i,j),0.0)/g);
    plk_n(i,j) = (td>0.0 && d->plk(i,j)>0.0) ? d->plk(i,j)*exp(-p->dt/td) : 0.0;
    plh_n(i,j) = plk_n(i,j)>1.0e-8 ? d->plh(i,j) : 0.0;
    
    if(plk_n(i,j)<=1.0e-8)
    plk_n(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,brkage_n,1);
    pgc->gcsl_start4(p,brkk_n,1);
    pgc->gcsl_start4(p,brkr_n,1);
    pgc->gcsl_start4(p,brkh_n,1);
    
    count_n=p->count;
    }
    return;
    }
    
    if(part==1)
    {
    const bool seeded = (p->A550==4 || p->A550==5);
    kn_count=0;
    pl_count=0;
    pl_xsum=0.0;
    on_count=0;
    on_xsum=0.0;
    
    if(p->A542>0.0)
    onset_filter(p,d,pgc,WL);
    
    SLICELOOP4
    {
    brkage(i,j) = 0.0;
    brkk(i,j) = 0.0;
    brkr(i,j) = 0.0;
    brkh(i,j) = 0.0;
    brkflag(i,j) = 0;
    d->vb(i,j) = 0.0;
    d->plk(i,j) = p->wet[IJ]==1 ? plk_n(i,j) : 0.0;
    d->plh(i,j) = p->wet[IJ]==1 ? plh_n(i,j) : 0.0;
    
        if(p->wet[IJ]==1 && brkzone(i,j)==0 && p->count>1)
        {
        const double h = MAX(WL(i,j),0.0);
        const double et = (p->A542>0.0) ? etf(i,j) : d->detadt(i,j);
        
        // the largest age of the column and its neighbours, and the wavenumber and initial threshold of that event (A 550 4)
        double age = brkage_n(i,j);
        double kw = brkk_n(i,j);
        double r0 = brkr_n(i,j);
        double hb = brkh_n(i,j);
        
        auto take = [&](int ii, int jj) { if(brkage_n(ii,jj)>age) { age = brkage_n(ii,jj); kw = brkk_n(ii,jj); r0 = brkr_n(ii,jj); hb = brkh_n(ii,jj); } };
        take(i-1,j);
        take(i+1,j);
        
        if(p->j_dir==1)
        {
        take(i,j-1);
        take(i,j+1);
        }
        
        // A 550 4: a seeded column starts or continues an event with the wavenumber of the breaking crest. Columns
        // without seed start events as in A 550 3 (deta/dt > A 554 c): with the seed as the only onset, the noise of
        // the crest speed C (B = u_s/C) moved the onset from wave to wave by up to 2.6 m in Ting & Kirby, and the wave
        // heights around the break point oscillated from gauge to gauge
        bool seedc = false;

        if(seeded && brkseed(i,j)==1)
        {
            if(age<=0.0 || kw<=0.0)
            kw = brkkseed(i,j);

            seedc = true;
        }
        
        // length of the wave: h (A 550 3), tanh(k h)/k (A 550 4, deep water 1/k); c = sqrt(g L_k)
        const double Lk = (seeded && kw>0.0) ? tanh(kw*h)/kw : h;
        const double c = sqrt(g*Lk);
        const double tstar = p->A503*sqrt(Lk/g);
        
        // initial threshold eta_t^I / c: A 554; a new event of a seed (A 550 4) starts at the deta/dt of the column if
        // that is larger, so that B starts from 0 and the viscosity grows continuously
        if(age<=0.0 || r0<=0.0)
        r0 = (seedc && c>0.0) ? MAX(p->A554, et/c) : p->A554;
        
        const double etI = r0*c;
        const double etF = p->A501*c;
        
        const double ets = (age>0.0 && tstar>0.0) ? etI + MIN(age/tstar,1.0)*(etF - etI) : etI;
        
        // S 137: the breaker height of a new event from its crest, H_b = eta_c/0.75 (crest of a breaker about 3/4 of
        // its height above the still water level); carried with the front as the age
        // no plunge of an event that starts in the swash or on a shoal shallower than H_b/2 (still-water depth WL - eta):
        // a wetting front has deta/dt > 0 and counts as breaking, but it is no plunging breaker (hb < 0 marks it)
        if(age<=0.0 || hb==0.0)
        {
        hb = MAX(d->eta(i,j),0.0)/0.75;
        
        if(WL(i,j) - d->eta(i,j) < 0.5*hb)
        hb = -1.0;
        }
        
            if(et>ets && ets>0.0)
            {
            const double B = (et>=2.0*ets) ? 1.0 : et/ets - 1.0;
            
            d->vb(i,j) = B*p->A504*p->A504*Lk*et;
            brkage(i,j) = age + p->dt;
            brkk(i,j) = kw;
            brkr(i,j) = r0;
            brkh(i,j) = hb;
            brkflag(i,j) = 1;
            ++kn_count;
            
            if(age<=0.0)
            {
            ++on_count;
            on_xsum += p->XP[IP];
            }
            
            plunge(p,d,i,j,age + p->dt,hb,kw,c);
            }
            
            // A 550 4: the event of a seeded column ages also below the threshold, so the threshold relaxes from
            // eta_t^I at the onset of the crest and the viscosity grows continuously over T*
            else
            if(seedc)
            {
            brkage(i,j) = age + p->dt;
            brkk(i,j) = kw;
            brkr(i,j) = r0;
            brkh(i,j) = hb;
            }
        }
    }
    
    // spread over the band of A 502 h (A 550 4: A 502 L_k of the event): bd is the largest width in cells
    kn_band=0;
    
    SLICELOOP4
    {
    double dxl = p->DXN[IP];
    
    if(p->j_dir==1)
    dxl = MIN(dxl,p->DYN[JP]);
    
    bd(i,j) = (p->wet[IJ]==1 && brkzone(i,j)==0) ? int(ceil(p->A502*MAX(WL(i,j),0.0)/dxl - 1.0e-9)) : 0;
    kn_band = MAX(kn_band,bd(i,j));
    }
    
    SLICELOOP4
    d->vbk(i,j) = brkk(i,j);
    return;
    }
    
    if(part==2)
    {
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4(p,d->vbk,1);
    return;
    }
    
    if(part==3)
    {
    const bool seeded = (p->A550==4 || p->A550==5);
    
    SLICELOOP4
    {
    brkvb(i,j) = d->vb(i,j);
    bt_vbs(i,j) = d->vbk(i,j);
    
        if(qn<bd(i,j))
        {
        // A 550 4: the band of an event is A 502 L_k wide, L_k = tanh(k h)/k with the wavenumber k of the event
        // carried with vb (<= h; with A 502 h in deep water the band reached over the next crest and trough)
        const double hb = MAX(WL(i,j),0.0);
        const double dxb = (p->j_dir==1) ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];
        auto inband = [&](double kv) { return !(seeded && kv>0.0) || qn < int(ceil(p->A502*tanh(kv*hb)/kv/dxb - 1.0e-9)); };
        auto take = [&](int ii, int jj) { if(d->vb(ii,jj)>brkvb(i,j) && inband(d->vbk(ii,jj))) { brkvb(i,j) = d->vb(ii,jj); bt_vbs(i,j) = d->vbk(ii,jj); } };
        take(i-1,j);
        take(i+1,j);
        
        if(p->j_dir==1)
        {
        take(i,j-1);
        take(i,j+1);
        }
        }
    }
    
    SLICELOOP4
    {
    d->vb(i,j) = brkvb(i,j);
    d->vbk(i,j) = bt_vbs(i,j);
    }
    return;
    }
    
    SLICELOOP4
    d->breaking(i,j) = brkflag(i,j);
    
    swash_taper(p,d,WL);
    
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4(p,d->vbk,1);
    pgc->gcsl_start4int(p,d->breaking,50);
    
    if(p->S137>0.0)
    {
    pgc->gcsl_start4(p,d->plk,1);
    pgc->gcsl_start4(p,d->plh,1);
    }
}

/*--------------------------------------------------------------------
A 542: onset of breaking independent of the grid (A 550 3, 4, 5). The test deta/dt > eta_t* (A 554 sqrt(g h)) sees
the front of a captured bore, which is 4-5 cells wide on every grid: deta/dt at the front grows as 1/dx, and the
onset moved with the grid (Baldock et al. 2011: 4 m between dx 0.1 and 0.05 m; Roelvink 1989 at dx 0.1 m never
reached A 554 0.6). The integral of deta/dt across a front, c H, does not depend on dx: deta/dt filtered over a
physical length l = A 542 h (Gaussian of the same variance as a box of width l, as explicit diffusion passes) is the
same on every grid finer than l. The onset test and the viscosity of breaking_kennedy use the filtered value.
Conservative passes (the dry columns and the relaxation zones are closed), sigma^2 = l^2/12.
--------------------------------------------------------------------*/
void nhflow_breaking::onset_filter(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    double r2max = 0.0;
    
    SLICELOOP4
    {
    etf(i,j) = (p->wet[IJ]==1 && brkzone(i,j)==0) ? d->detadt(i,j) : 0.0;
    
        if(p->wet[IJ]==1)
        {
        const double l = p->A542*MAX(WL(i,j),0.0);
        const double dx = (p->j_dir==1) ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];
        r2max = MAX(r2max, l*l/(dx*dx));
        }
    }
    
    r2max = pgc->globalmax(r2max);
    
    // passes: variance 2 a dx^2 per pass and direction, a <= 0.25 (stable in 2D), total l^2/12
    const int np = MAX(1, int(ceil(r2max/(12.0*2.0*0.25))));
    
    pgc->gcsl_start4(p,etf,1);
    
    for(int q=0; q<np; ++q)
    {
        SLICELOOP4
        {
        etf2(i,j) = etf(i,j);
        
            if(p->wet[IJ]==1 && brkzone(i,j)==0)
            {
            const double l = p->A542*MAX(WL(i,j),0.0);
            
            auto coef = [&](int ii, int jj, double dx)
            {
                if(p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]!=1 || brkzone(ii,jj)!=0)
                return 0.0;
                
                const double ln = 0.5*(l + p->A542*MAX(WL(ii,jj),0.0));
                
                return MIN(0.25, ln*ln/(24.0*double(np)*dx*dx));
            };
            
            double sum = 0.0;
            
            const double ce = coef(i+1,j,p->DXP[IP]);
            const double cw = coef(i-1,j,p->DXP[IM1]);
            sum += ce*(etf(i+1,j)-etf(i,j)) + cw*(etf(i-1,j)-etf(i,j));
            
                if(p->j_dir==1)
                {
                const double cn = coef(i,j+1,p->DYP[JP]);
                const double cs = coef(i,j-1,p->DYP[JM1]);
                sum += cn*(etf(i,j+1)-etf(i,j)) + cs*(etf(i,j-1)-etf(i,j));
                }
            
            etf2(i,j) = etf(i,j) + sum;
            }
        }
        
        SLICELOOP4
        etf(i,j) = etf2(i,j);
        
        pgc->gcsl_start4(p,etf,1);
    }
}

/*--------------------------------------------------------------------
S 137: impact of the plunging jet (A 550 3, 4, 5). A single-valued surface turns the plunging jet into a bore, and
the breaking viscosity is spread over a band behind the front, so nothing marks the point where the jet hits the
water and drives the turbulence down to the bed (the trough under a plunging breaker, Baldock et al. 2011). The
jet lands a plunge distance l_p = S 138 H_b shoreward of the onset (Galvin 1969, Smith & Kraus 1991: 1.5 - 4 H_b).
The front travels with the event of the breaking model: its age t_b is the time since the onset, so the front is
at the impact point when t_b = l_p/c, and it crosses a footprint of width H_b in H_b/c. In that window the column
carries the impact turbulence at the surface,
    k_p = S 137 c^2        (the jet falls with about the crest speed c; S 137 its share as turbulent energy),
decaying with the depth below the surface as exp(-d/(S 139 H_b)), and H_b of the event (d->plk, d->plh). With
S 142 > 0 the impact turbulence stays in the column after the window and decays over S 142 sqrt(H_b/g). The bed
shear (bedshear_nhflow, u_eff^2 = u^2 + 1.45 k_b) and the suspended-load diffusivity (nhflow_suspended_IM1,
sqrt(k) 0.1 H_b) take it. H_b = eta_c/0.75 from the crest at the onset; no plunge for an event that starts where the
still-water depth is below H_b/2 (swash, wetting fronts), and no impact where the still-water depth is below H_b/4.
Breaker type: with the wavenumber k of the crest (A 550 4, 5) only plunging breakers, Iribarren number
xi_b = tan(beta)/sqrt(H_b/L_0) >= 0.4 (Battjes 1974), L_0 = 2 pi/(k tanh(k h)); A 550 3 has no k: every event.
--------------------------------------------------------------------*/
void nhflow_breaking::plunge(lexer *p, fdm_nhf *d, int ii, int jj, double age, double hb, double kw, double c)
{
    if(p->S137<=0.0 || hb<=1.0e-4 || c<=0.0)
    return;   // hb < 0: an event of the swash (no plunge)
    
    const double tp = p->S138*hb/c;
    const double tw = hb/c;
    
    if(age<tp || age>=tp+tw)
    return;
    
    const int i0=i, j0=j;
    i=ii; j=jj;
    
    // the jet lands in water: still-water depth at the impact at least H_b/4 (not in the swash)
    bool plunging = d->WL(i,j) - d->eta(i,j) >= 0.25*hb;
    
    if(plunging && kw>0.0)
    {
    const double h = MAX(d->WL(i,j),1.0e-6);
    double sx = (d->bed(i+1,j)-d->bed(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
    double sy = p->j_dir==1 ? (d->bed(i,j+1)-d->bed(i,j-1))/(p->DYP[JP]+p->DYP[JM1]) : 0.0;
    const double L0 = 6.283185307179586/(kw*tanh(kw*h));
    const double xi = sqrt(sx*sx + sy*sy)/sqrt(hb/L0);
    
    plunging = xi>=0.4;
    }
    
    if(plunging && p->S137*c*c>d->plk(i,j))
    {
    d->plk(i,j) = p->S137*c*c;
    d->plh(i,j) = hb;
    ++pl_count;
    pl_xsum += p->XP[IP];
    }
    
    i=i0; j=j0;
}

// A 553: no breaking viscosity in thin water, h < A 553, linear to the full value at 2 A 553 (as S 132 for the
// sediment). The swash and the run-up tongue have a rising level (deta/dt > 0) and a threshold that goes to zero
// with the depth (A 554 sqrt(g h)), so every wetting front counts as breaking; A 550 1 then puts the full A 557
// into a layer of a few mm. The flags (output, breaking count) are kept.
void nhflow_breaking::swash_taper(lexer *p, fdm_nhf *d, slice &WL)
{
    if(p->A553<=0.0)
    return;
    
    SLICELOOP4
    d->vb(i,j) *= MAX(0.0, MIN(1.0, MAX(WL(i,j),0.0)/p->A553 - 1.0));
}

// the relaxation zones of the wave generation and the numerical beach (ioflow::relax_boxes): deta/dt of the
// continuity is not the change of the relaxed water level there
void nhflow_breaking::breaking_zone_mask(lexer *p, const std::function<void(const std::function<void(int,int)>&)> &cells)
{
    SLICELOOP4
    brkzone(i,j) = 0;

    cells([&](int ii, int jj) { brkzone(ii,jj) = 1; });

    brkzone_built = true;
}

void nhflow_breaking::breaking_zones(lexer *p, ioflow *pflow)
{
    if(brkzone_built || p->A550<2 || pflow==nullptr)
    return;
    
    // boxes xs, xe, ys, ye
    std::vector<double> b;
    pflow->relax_boxes(p,b);
    
    SLICELOOP4
    {
    const double x = (p->geo_curv==1) ? p->XC2D[(i-p->imin)*p->jmax + (j-p->jmin)] : p->XP[IP];
    const double y = (p->geo_curv==1) ? p->YC2D[(i-p->imin)*p->jmax + (j-p->jmin)] : p->YP[JP];
    
    brkzone(i,j) = 0;
    
        for(size_t q=0; q+3<b.size(); q+=4)
        if(x>=b[q] && x<=b[q+1] && (p->j_dir==0 || (y>=b[q+2] && y<=b[q+3])))
        brkzone(i,j) = 1;
    }
    
    brkzone_built = true;
}

void nhflow_breaking::breaking_baquet(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    // mesh refinement (nhflow_amr, brk_defer): the three parts run on all grids in turn, the
    // patches exchange the detection and the flags of their sibling cells in between
    if(brk_defer)
    {
    brk_eta = &eta;
    brk_eta_n = &eta_n;
    brk_WL = &WL;
    brk_alpha = alpha;
    return;
    }

    breaking_detect(p,d,pgc,eta,eta_n,WL,alpha);
    breaking_flags(p,d,pgc);
    breaking_vb(p,d,pgc,WL);
}

void nhflow_breaking::breaking_part(lexer *p, fdm_nhf *d, ghostcell *pgc, int part)
{
    if(part==0)
    breaking_detect(p,d,pgc,*brk_eta,*brk_eta_n,*brk_WL,brk_alpha);
    if(part==1)
    breaking_flags(p,d,pgc);
    if(part==2)
    breaking_vb(p,d,pgc,*brk_WL);
}

// detection: steepness (bx, by: 10 / 20) and depth induced (bd: 2)
void nhflow_breaking::breaking_detect(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    if(p->count>count_n)
    {
    SLICELOOP4
    {
    brkflag(i,j)=0;
    d->breaking(i,j)=0;
    }
    
    count_n=p->count;
    }
    
    SLICELOOP4
    {
    bx(i,j)=0;
    by(i,j)=0;
    bd(i,j)=0;
    }
    
    pgc->gcsl_start4int(p,brkflag,50);
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    pgc->gcsl_start4int(p,bd,50);
    
    // steepness induced breaking
    if((p->A551==2 || p->A551==3) && p->count>1)
    {
    SLICELOOP4
    {
            // x
            if(d->Ex(i,j)   < -p->A555)
            {
                bx(i,j) = 10;
            }
            
            if(d->Ex(i,j)   > p->A555)
            {
                bx(i,j) = 20;

            }
            
            // y
            if(p->j_dir==1)
            if( d->Ey(i,j)   < -p->A555)
            {
                by(i,j) = 10;
            }
            
            if(p->j_dir==1)
            if( d->Ey(i,j)   > p->A555)
            {
                by(i,j) = 20;
            }
            
    }
    
    // mesh refinement: no detection in the filled cells around a patch (their surface comes
    // from the coarser grid and their outer cells see one-sided fluxes); their breaking comes
    // from the source grid below
    if(p->amrvb!=nullptr)
    SLICELOOP4
    if(p->amrvb[IJ]>=0.0)
    {
    bx(i,j)=0;
    by(i,j)=0;
    }

    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    }
    
    // depth induced breaking
    if((p->A551==1 || p->A551==3) && p->count>1)
    {
    // deta/dt of the stage from the continuity (nhflow_fsf rk2/rk3 steps); (eta - eta_n)/(alpha dt) is right
    // only in the first stage (eta_n is the level of the previous stage)
    const double g = fabs(p->W22);
    
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(d->detadt(i,j) > p->A554*sqrt(g*MAX(WL(i,j),0.0)))
        bd(i,j)=2;
    }
    
    if(p->amrvb!=nullptr)
    SLICELOOP4
    if(p->amrvb[IJ]>=0.0)
    bd(i,j)=0;

    pgc->gcsl_start4int(p,bd,50);
    }
}

// the breaking flags from the detection
void nhflow_breaking::breaking_flags(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if((p->A551==2 || p->A551==3) && p->count>1)
    {
    SLICELOOP4
    if(bx(i,j)>0 || by(i,j)>0)
    {
    brkflag(i,j)=1;
    }
    
    SLICELOOP4
    {
    // x
    if(bx(i,j)==10 || bx(i-1,j)==10 || bx(i+1,j)==10 || bx(i+2,j)==10)
    brkflag(i,j)=1;
    
    if(bx(i,j)==20 || bx(i-1,j)==20 || bx(i-2,j)==20 || bx(i+1,j)==20)
    brkflag(i,j)=1;
    
    // y (by is zero with one cell in y)
    if(p->j_dir==1)
    if(by(i,j)==10 || by(i,j-1)==10 || by(i,j+1)==10 || by(i,j+2)==10)
    brkflag(i,j)=1;
    
    if(p->j_dir==1)
    if(by(i,j)==20 || by(i,j-1)==20 || by(i,j-2)==20 || by(i,j+1)==20)
    brkflag(i,j)=1;
        
    }
    }
    
    if((p->A551==1 || p->A551==3) && p->count>1)
    SLICELOOP4
    {
    // x
    if(bd(i,j)==2 || bd(i-1,j)==2 || bd(i-2,j)==2 || bd(i+1,j)==2 || bd(i+2,j)==2)
    brkflag(i,j)=2;
    
    // y
    if(p->j_dir==1)
    if(bd(i,j)==2 || bd(i,j-1)==2 || bd(i,j-2)==2 || bd(i,j+1)==2 || bd(i,j+2)==2)
    brkflag(i,j)=2;    
    }
    
    // mesh refinement: a filled cell breaks where its source cell carries the full breaking
    // viscosity
    if(p->amrvb!=nullptr)
    SLICELOOP4
    if(p->amrvb[IJ]>=p->A557 && brkflag(i,j)==0)
    brkflag(i,j)=1;
}

// the breaking viscosity from the flags
void nhflow_breaking::breaking_vb(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    // ------------------------
    // fill breaking viscosity
    // ------------------------
    
    SLICELOOP4
    d->vb(i,j) = 0.0;
        
    if(p->j_dir==0)
    SLICELOOP4
    {   
        // the cell and its direct neighbours A 557, the next ones A 557/2 (before, the second test
        // overwrote the direct neighbours with A 557/2)
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0)
        d->vb(i,j) = p->A557;
            
        else
        if(brkflag(i-2,j)>0 || brkflag(i+2,j)>0)
        d->vb(i,j) = 0.5*p->A557;
    }

    if(p->j_dir==1)
    SLICELOOP4
    {   
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0 || brkflag(i,j-1)>0 || brkflag(i,j+1)>0)
        d->vb(i,j) = p->A557;
            
        else
        if(( brkflag(i-1,j-1)>0 || brkflag(i-1,j+1)>0 || brkflag(i+1,j-1)>0 || brkflag(i+1,j+1)>0
        || brkflag(i-2,j)>0 || brkflag(i+2,j)>0 || brkflag(i,j-2)>0 || brkflag(i,j+2)>0))
        d->vb(i,j) = 0.5*p->A557;
    }
    
        
        // additional breaking filter
        // shallow
        if(p->A552==1)
        SLICELOOP4
        if(brkflag(i,j)==2)
        {
         filter(p,d,pgc,WL);
        }   
        
        // deep
        if(p->A552==2)
        SLICELOOP4
        if(brkflag(i,j)==1)
        {
         filter(p,d,pgc,WL);
        }   
        
        // all
        if(p->A552==3)
        SLICELOOP4
        if(brkflag(i,j)>=1)
        {
         filter(p,d,pgc,WL);
        }   
        
    // the flags for the output (VTP 'breaking') and the breaking count of the log
    SLICELOOP4
    d->breaking(i,j) = brkflag(i,j);
    
    swash_taper(p,d,WL);
    
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4int(p,d->breaking,50);
    
    SLICELOOP4
    d->breaklog(i,j)=0;
    
    // breaklog
    int count=0;
    
    SLICELOOP4
    if(d->breaking(i,j)>0)
    {
    d->breaklog(i,j)=1;
    ++count;
    }

    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}
