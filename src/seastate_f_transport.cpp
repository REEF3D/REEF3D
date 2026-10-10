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
#include"seastate_param.h"
#include"seastate_exchange.h"
#include"seastate_implicit.h"
#include"seastate_swan_spc.h"
#include"seastate_source.h"
#include"seastate_amr.h"
#include"sliceint4.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<sys/stat.h>

void seastate_f::boundary(lexer *p, ghostcell *pgc)
{
    side[0] = p->A712_xm;
    side[1] = p->A712_xp;
    side[2] = p->A712_ym;
    side[3] = p->A712_yp;

    Nb.clear();

    // surfbeat: boundary time series from the input spectrum (seastate_f_surfbeat.cpp)
    if(p->A770==1)
    {
    surfbeat_ini(p,pgc);
    return;
    }

    // time series of spectra at many locations (seastate_f_forcing.cpp)
    if(p->A711==3)
    {
    boundary_series(p,pgc);
    return;
    }

    boundary_spectrum(p,pgc,*e->grid,Nb);
}

void seastate_f::boundary_spectrum(lexer *p, ghostcell *pgc, const seastate_grid &g, std::vector<float> &N)
{
    if(p->A711==1)
    parametric_spectrum(p,pgc,g,N,"A 711 1");

    // SWAN 2D spectrum file, first location and time
    if(p->A711==2)
    {
    seastate_swan_spc spc;
    std::string err;

        if(!spc.read("seastate-boundary.spc",err))
        {
            if(p->mpirank==0)
            cout<<endl<<"SEASTATE input error  --  A 711 2: "<<err<<endl<<endl;

            pgc->final(true);
        }

    spc.to_grid(g,N);

    seastate_param sp;
    sp.compute(g,N.data());

        if(p->mpirank==0)
        cout<<"SEASTATE boundary spectrum (A 711 2): seastate-boundary.spc, "<<spc.f.size()<<" frequencies "<<spc.f.front()<<" - "<<spc.f.back()
            <<" Hz, "<<spc.dir.size()<<" directions; on the spectral grid: Hs "<<sp.Hs<<" m, Tp "<<sp.Tp<<" s, direction "<<sp.dir<<" deg"<<endl<<endl;
    }
}

// stationary convergence test of a cell: the change d of Hs in the last iteration (A 799 0), or (A 799 1)
// the estimated distance of Hs to the solution d/(1-rho), rho = d/d_previous the contraction of the
// iteration in the cell (at most 0.99; an alternating change counts with d); relative to Hs
static double seastate_distance(int test, double hs, double hs0, double &dprev, double floor)
{
    const double d = hs-hs0;
    double c = std::fabs(d);

    if(test==1 && d*dprev>0.0)
    c /= 1.0-std::min(std::fabs(d)/std::fabs(dprev),0.99);

    dprev = d;

    return c/std::max(hs,floor);
}

void seastate_f::transport(lexer *p, ghostcell *pgc)
{
    if(pamr!=nullptr)
    {
    transport_amr(p,pgc);
    return;
    }

    const bool stationary = (p->A700==2);

    // stationary: 1/dt = 0; with the deep-water physics pseudo time step A 706 (the wind-sea
    // source terms are lagged by one iteration, the pseudo time step damps the iteration)
    const double rdt = stationary ? ((psrc!=nullptr && psrc->param().komen) ? 1.0/p->A706 : 0.0) : 1.0/dtw;
    const bool refraction = (p->A713==1);
    const bool fshift = (p->A714==1 && sb==nullptr);      // surfbeat: one frequency, no frequency shift

    iter_done = 0;
    conv = 0.0;

    if(!stationary)
    {
        if(N0!=nullptr)
        N0->copy_from(*e->N);

        for(int it=0; it<iter_max; ++it)
        {
        diffraction(p,pgc);
        obstacles(p,pgc);
        psolv->iterate(p,pgc,e,pex,N0,rdt,Nb,side,refraction,fshift);
        ++iter_done;
        }

        return;
    }

    // stationary: iterate until the largest relative change of Hs is below A 708
    const seastate_grid &g = *e->grid;
    seastate_param sp;

    hs_old.assign(size_t(p->knox)*p->knoy,0.0);
    std::vector<double> dprev(hs_old.size(),0.0);

    SLICELOOP4
    if(e->wet(i,j)==1)
    {
    sp.compute(g,e->N->spec(i,j));
    hs_old[size_t(i)*p->knoy+j] = sp.Hs;
    }

    for(int it=0; it<iter_max; ++it)
    {
    // Phase 6: diffraction parameter and obstacle transmission from the latest spectra (A 718, A 722)
    diffraction(p,pgc);
    obstacles(p,pgc);
    psolv->iterate(p,pgc,e,pex,nullptr,rdt,Nb,side,refraction,fshift);
    ++iter_done;

    double hmax=0.0;
    std::vector<double> hs(size_t(p->knox)*p->knoy,0.0);

        SLICELOOP4
        if(e->wet(i,j)==1)
        {
        sp.compute(g,e->N->spec(i,j));
        hs[size_t(i)*p->knoy+j] = sp.Hs;
        hmax = std::max(hmax,sp.Hs);
        }

    hmax = pgc->globalmax(hmax);

    double cmax=0.0;
    const double floor = std::max(0.01*hmax,1.0e-6);
    double nok=0.0, nall=0.0;

        SLICELOOP4
        if(e->wet(i,j)==1)
        {
        const size_t n = size_t(i)*p->knoy+j;
        const double c = seastate_distance(p->A799,hs[n],hs_old[n],dprev[n],floor);
        cmax = std::max(cmax,c);
        nall += 1.0;
        if(c<p->A708)
        nok += 1.0;
        }

    conv = pgc->globalmax(cmax);
    hs_old.swap(hs);

    nok = pgc->globalsum(nok);
    nall = pgc->globalsum(nall);
    const double percent = nall>0.0 ? 100.0*nok/nall : 100.0;

        // convergence history (rank 0): REEF3D_SEASTATE_Log/REEF3D_SEASTATE_convergence.dat
        if(p->mpirank==0)
        {
            if(!convlog.is_open())
            {
            mkdir("./REEF3D_SEASTATE_Log",0777);
            convlog.open("./REEF3D_SEASTATE_Log/REEF3D_SEASTATE_convergence.dat");
            convlog<<"# count \t iteration \t max. relative change of Hs \t percentage of the cells within A 708"<<endl;
            }
        convlog<<p->count<<" \t "<<iter_done<<" \t "<<setprecision(6)<<conv<<" \t "<<percent<<endl;
        }

        if(conv<p->A708)
        break;

        // SWAN-type criterion (A 709 < 100): the percentage of cells within A 708
        if(p->A709<100.0 && percent>=p->A709)
        break;
    }
}

// mesh refinement (G 1): composite iterations over all grids (seastate_amr); stationary runs iterate
// until the largest relative change of Hs on the leaf cells of all grids is below A 708, or (A 709)
// the given percentage of the area of the leaf cells meets A 708 (cells weighted with their area, so
// that refining the shallow and coastal cells does not change the weight of the criterion)
void seastate_f::transport_amr(lexer *p, ghostcell *pgc)
{
    const bool stationary = (p->A700==2);
    const double rdt = stationary ? ((psrc!=nullptr && psrc->param().komen) ? 1.0/p->A706 : 0.0) : 1.0/dtw;
    const bool refraction = (p->A713==1);
    const bool fshift = (p->A714==1);

    iter_done = 0;
    conv = 0.0;

    if(!stationary)
    {
        if(N0!=nullptr)
        N0->copy_from(*e->N);

        pamr->step_begin();

        for(int it=0; it<iter_max; ++it)
        {
        diffraction(p,pgc);
        obstacles(p,pgc);
        pamr->iterate(p,pgc,N0,rdt,Nb,side,refraction,fshift);
        ++iter_done;

        // FAS coarse-grid correction (A 758) within the step: every A 758 n iterations of the step, not after the
        // last one (the prolonged correction needs an iteration of the fine grids after it, else its error enters
        // the next step: with one iteration per step and a correction after it, the wind sea grew 18 % too fast)
        if(p->A758_n>0 && iter_done%p->A758_n==0 && it<iter_max-1)
        pamr->fas(p,pgc,N0,rdt,Nb,side,refraction,fshift,p->A758_m);
        }

        return;
    }

    const seastate_grid &g = *e->grid;
    seastate_param sp;

    // Hs of the leaf cells: level 0 without the covered cells, then the patches; area weights
    std::vector<double> wt;
    auto leaf_hs = [&](std::vector<double> &hs, double &hmax)
    {
        hs.clear();
        wt.clear();
        hmax = 0.0;

        SLICELOOP4
        if(e->wet(i,j)==1 && pamr->level0_covered(i,j)==0)
        {
        sp.compute(g,e->N->spec(i,j));
        hs.push_back(sp.Hs);
        wt.push_back(1.0);
        hmax = std::max(hmax,sp.Hs);
        }

        double vmin = 0.0;
        pamr->parameters(&hs,hmax,vmin,&wt);
    };

    std::vector<double> hs, hs0;
    double hmax;
    leaf_hs(hs0,hmax);
    std::vector<double> dprev(hs0.size(),0.0);

    for(int it=0; it<iter_max; ++it)
    {
    // Phase 6: diffraction, obstacles and coasts on all grids (A 718, A 722 - A 726)
    diffraction(p,pgc);
    obstacles(p,pgc);
    pamr->iterate(p,pgc,nullptr,rdt,Nb,side,refraction,fshift);
    ++iter_done;

    // FAS coarse-grid correction (A 758): level 0 as the coarse grid of the patches
    if(p->A758_n>0 && iter_done%p->A758_n==0)
    {
    pamr->fas(p,pgc,nullptr,rdt,Nb,side,refraction,fshift,p->A758_m);
    if(p->mpirank==0)
    cout<<"SEASTATE FAS: iteration "<<iter_done<<", "<<p->A758_m<<" coarse iterations, largest relative change of the coarse cells "<<pamr->fas_last()<<endl;
    }

    leaf_hs(hs,hmax);
    hmax = pgc->globalmax(hmax);

    double cmax=0.0;
    const double floor = std::max(0.01*hmax,1.0e-6);

    double nok=0.0, nall=0.0;

        for(size_t n=0; n<hs.size() && n<hs0.size(); ++n)
        {
        if(dprev.size()<=n)
        dprev.resize(n+1,0.0);

        const double c = seastate_distance(p->A799,hs[n],hs0[n],dprev[n],floor);
        cmax = std::max(cmax,c);
        nall += wt[n];
        if(c<p->A708)
        nok += wt[n];
        }

    conv = pgc->globalmax(cmax);
    hs0.swap(hs);
    nok = pgc->globalsum(nok);
    nall = pgc->globalsum(nall);

    pamr->convergence(p,iter_done,conv,nall>0.0 ? 100.0*nok/nall : 100.0);

        if(conv<p->A708)
        break;

        if(p->A709<100.0 && nall>0.0 && 100.0*nok/nall>=p->A709)
        break;
    }
}
