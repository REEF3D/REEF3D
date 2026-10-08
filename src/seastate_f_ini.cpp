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
#include"seastate_vtp.h"
#include"seastate_exchange.h"
#include"seastate_implicit.h"
#include"seastate_source.h"
#include"seastate_roller.h"
#include"seastate_amr.h"
#include"seastate_bathy.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<sys/stat.h>
#include<sys/types.h>

void seastate_f::ini(lexer *p, ghostcell *pgc)
{
    ini_common(p,pgc,false);
}

void seastate_f::ini_coupled(lexer *p, ghostcell *pgc)
{
    ini_common(p,pgc,true);
}

void seastate_f::ini_common(lexer *p, ghostcell *pgc, bool coupled_)
{
    coupled = coupled_;

    check_keys(p,pgc);

    dtw = p->A706;

    if(!coupled)
    {
    p->count=0;
    p->printcount=0;
    p->dt=p->A706;

    // 2D points and polygons for the output (as SFLOW)
    int count=0;
    p->pointnum2D=0;
    p->cellnum2D=0;
    p->polygon_sum=0;

    TPSLICELOOP
    {
    ++count;
    ++p->pointnum2D;
    e->nodeval(i,j)=count;
    }

    SLICEBASELOOP
    ++p->polygon_sum;

    p->polygon_sum*=2;

    SLICELOOP4
    ++p->cellnum2D;

    p->cellnumtot2D=pgc->globalisum(p->cellnum2D);
    }
    else
    {
    // the host has set the point and polygon counts (same 2D grid, same numbering as SFLOW)
    int count=0;

    TPSLICELOOP
    {
    ++count;
    e->nodeval(i,j)=count;
    }
    }

    environment(p,pgc);
    storage(p,pgc);
    kinematics(p,pgc,0.0);
    initial(p,pgc);

    if(p->A711==3 || p->A730==2)
    forcing_ini(p,pgc);

    boundary(p,pgc);

    // transport
    // nonstationary: 2 iterations on several ranks, so that the lagged halo values of the first
    // sweep are corrected (energy conservation across rank borders); surfbeat: always 2 (the
    // second-order correction is lagged, results independent of the decomposition)
    iter_max = p->A707>0 ? p->A707 : (p->A700==2 ? 50 : ((p->M10>1 || p->A770==1) ? 2 : 1));

    // surfbeat with second-order advection (A 775 2): 2 halo layers
    const bool second = (p->A770==1 && p->A775==2);

    pex   = new seastate_exchange(p,e->grid->nbin,(second || p->A796==2) ? 2 : 1);
    psolv = new seastate_implicit(p,e);
    psolv->sparsity(p->A795);
    psolv->geographic_order(p->A796);
    psolv->threads(p->A798);
    if(p->A700==2)
    psolv->source_iterations(p->A738,p->A739);

    if(sb!=nullptr)
    psolv->boundary_rows(&Nbx,&Nbx0);

    if(bser!=nullptr)
    {
    const std::vector<float> *s[4] = {&Nside[0],&Nside[1],&Nside[2],&Nside[3]};
    psolv->boundary_sides(s);
    }

    if(wser!=nullptr)
    psolv->wind_field(wU10,wdir);

    if(bser!=nullptr || wser!=nullptr)
    forcing_update(p,pgc,p->simtime);

    psolv->second_order(second);

    sources(p,pgc);

    // surfbeat roller (source: Roelvink breaking)
    if(p->A770==1 && p->A748==1 && p->A740==2)
    proll = new seastate_roller(p,e,p->A749);

    if(p->A700==1 && iter_max>1)
    {
    N0 = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,e->grid->nbin,p->A704);
    N0->build(e->wet.V);
    }

    pex->start(p,pgc,*e->N);

    // mesh refinement (G 1): static hierarchy of patches
    if(p->G1>0)
    {
    pamr = new seastate_amr(p,pgc);
    seastate_amr::level0 l0 = {e,psolv,pex,psrc,N0,bathy,wser,tref};
    pamr->ini(p,pgc,l0);

        if(wser!=nullptr)
        pamr->wind(tref+p->simtime);
    }

    parameters(p,pgc);

    pprint = new seastate_vtp(p,e,pgc,coupled);

    if(wser!=nullptr)
    {
    pprint->add_field("U10x",wUx);
    pprint->add_field("U10y",wUy);
    }

    if(proll!=nullptr)
    {
    pprint->add_field("roller",&proll->R);
    pprint->add_field("Dr",&proll->Dr);
    pprint->add_field("Dw",&proll->Dw);
    }

    log_ini(p);

    // initial state (coupled: printed by the host after it has added its fields)
    if(!coupled)
    {
    const int pc = p->printcount;
    pprint->start(p,e,pgc);
    if(pamr!=nullptr && p->printcount!=pc)
    pamr->print(p,pgc);
    }

    log_step(p);

    if(coupled && p->A760>0)
    handover(p,pgc);
}

void seastate_f::check_keys(lexer *p, ghostcell *pgc)
{
    const char *msg = nullptr;

    if(p->A700!=1 && p->A700!=2)
    msg = "A 700: mode must be 1 (nonstationary) or 2 (stationary)";
    else if(p->A704<1)
    msg = "A 704: the tile size must be at least 1";
    else if(p->A705<0.0)
    msg = "A 705: the minimum water depth must not be negative";
    else if(!(p->A706>0.0))
    msg = "A 706: the time step must be positive";
    else if(p->A707<0)
    msg = "A 707: the number of iterations must not be negative";
    else if(!(p->A708>0.0))
    msg = "A 708: the convergence criterion must be positive";
    else if(p->A710!=0 && p->A710!=1)
    msg = "A 710: initial spectrum must be 0 (zero) or 1 (parametric)";
    else if(p->A711<0 || p->A711>3)
    msg = "A 711: boundary spectrum must be 0 (none), 1 (parametric), 2 (SWAN spectrum file) or 3 (time series of spectra at many locations)";
    else if(p->A712_xm<0 || p->A712_xm>2 || p->A712_xp<0 || p->A712_xp>2 || p->A712_ym<0 || p->A712_ym>2 || p->A712_yp<0 || p->A712_yp>2)
    msg = "A 712: each side must be 0 (open), 1 (boundary spectrum) or 2 (zero gradient)";
    else if(p->A720!=0 && p->A720!=1)
    msg = "A 720: prescribed current must be 0 (none) or 1 (linear in x, A 721)";
    else if(p->A730<0 || p->A730>2)
    msg = "A 730: wind must be 0 (none), 1 (uniform, A 731) or 2 (field, seastate-wind.dat)";
    else if(p->A730==1 && !(p->A731_u10>=0.0))
    msg = "A 731: the wind speed must not be negative";
    else if(p->A730>=1 && p->A732!=1)
    msg = "A 730: wind input needs the deep-water physics A 732 1 (Komen)";
    else if(p->A732!=0 && p->A732!=1)
    msg = "A 732: deep-water physics must be 0 (off) or 1 (Komen)";
    else if(p->A733!=0 && p->A733!=1)
    msg = "A 733: quadruplets must be 0 (off) or 1 (DIA)";
    else if(!(p->A734>=0.0))
    msg = "A 734: the linear growth coefficient must not be negative";
    else if(!(p->A735>=0.0))
    msg = "A 735: the limiter coefficient must not be negative";
    else if(p->A740!=0 && p->A740!=1 && p->A740!=2)
    msg = "A 740: depth-induced breaking must be 0 (off), 1 (Battjes-Janssen) or 2 (Roelvink, surfbeat)";
    else if(p->A740==2 && p->A770!=1)
    msg = "A 740 2: Roelvink breaking is the breaking of the surfbeat model (A 770 1)";
    else if(p->A740==1 && p->A770==1)
    msg = "A 740 1: the surfbeat model (A 770 1) needs Roelvink breaking (A 740 2) or none";
    else if(p->A740>=1 && !(p->A741_alpha>0.0 && p->A741_gamma>0.0))
    msg = "A 741: the breaking coefficients alpha and gamma must be positive";
    else if(p->A742!=0 && p->A742!=1)
    msg = "A 742: bottom friction must be 0 (off) or 1 (JONSWAP)";
    else if(!(p->A743>=0.0))
    msg = "A 743: the friction coefficient must not be negative";
    else if(p->A744!=0 && p->A744!=1)
    msg = "A 744: triads must be 0 (off) or 1 (LTA)";
    else if(!(p->A745>=0.0))
    msg = "A 745: the triad coefficient must not be negative";
    else if(!(p->A736_cutfr>1.0) || !(p->A736_urcrit>0.0) || !(p->A736_urslim>=0.0))
    msg = "A 736: the triad parameters cutfr (> 1), urcrit (> 0) and urslim (>= 0) are out of range";
    else if(p->A737!=0 && p->A737!=1)
    msg = "A 737: the maximum energy must be 0 (off) or 1 (on)";
    else if(p->A750==1 && p->A10!=2 && p->A10!=5)
    msg = "A 750 1: the coupling needs SFLOW (A 10 2) or NHFLOW (A 10 5) as the host model";
    else if(p->A750==1 && p->A751!=1 && p->A751!=2)
    msg = "A 751: the wave forcing must be 1 (radiation stress) or 2 (vortex force)";
    else if(p->A750==1 && (p->A753<0 || p->A753>2))
    msg = "A 753: the feedback must be 0 (none), 1 (water level) or 2 (water level and currents)";
    else if(!(p->A752>=0.0))
    msg = "A 752: the ramp-up time must not be negative";
    else if(!(p->A761>0.0))
    msg = "A 761: the half-width of the handover sector must be positive";
    else if(p->A770!=0 && p->A770!=1)
    msg = "A 770: surfbeat must be 0 (off) or 1 (on)";
    else if(p->A770==1 && p->A700!=1)
    msg = "A 770 1: the surfbeat model is nonstationary (A 700 1)";
    else if(p->A770==1 && (p->A730!=0 || p->A732!=0 || p->A744!=0))
    msg = "A 770 1: wind input, whitecapping, quadruplets and triads need a frequency spectrum (A 730 0, A 732 0, A 744 0)";
    else if(p->A770==1 && p->A760>0)
    msg = "A 770 1: the handover (A 760) needs a frequency spectrum, not the wave groups";
    else if(p->A770==1 && p->A710!=0)
    msg = "A 770 1: the surfbeat model starts from rest (A 710 0)";
    else if(p->A780<0.0 || (p->A780>0.0 && p->A780<10000101.0))
    msg = "A 780: the start date-time must be YYYYMMDD.HHMMSS (or 0)";
    else if(!(p->A709>0.0) || p->A709>100.0)
    msg = "A 709: the percentage of converged cells must be in (0, 100]";
    else if(p->A790!=0 && p->A790!=1)
    msg = "A 790: the bathymetry raster must be 0 (off) or 1 (seastate-bathy.dat)";
    else if(p->A790==1 && coupled)
    msg = "A 790 1: the bathymetry raster is for stand-alone runs (A 10 7); a host model sets the bed";
    else if(p->G1>0 && coupled)
    msg = "G 1: mesh refinement of SEASTATE is for stand-alone runs (A 10 7), not with the coupling (A 750)";
    else if(p->G1>0 && p->A770==1)
    msg = "G 1: mesh refinement is not available with the surfbeat model (A 770 1)";
    else if(p->G1>0 && p->G40!=0)
    msg = "G 40: SEASTATE patches are cut at the rank boxes (G 40 0)";
    else if(p->G1>0 && (p->G4<4 || p->G4%2!=0))
    msg = "G 4: the tile size must be even and at least 4";
    else if(p->A791<0.0 || p->A792<0 || p->A792>3 || p->A793<0.0 || p->A794<0)
    msg = "A 791-794: the refinement criteria must not be negative (A 792 at most 3)";
    else if(p->A770==1 && p->A711!=1 && p->A711!=2)
    msg = "A 770 1: the surfbeat model needs a boundary spectrum (A 711 1 parametric or 2 SWAN file)";
    else if(p->A770==1 && p->A712_xm!=1)
    msg = "A 770 1: the wave groups enter through the x- side (A 712 1 ...)";
    else if(p->A770==1 && !(p->A772>0.0))
    msg = "A 772: the record length must be positive";
    else if(p->A770==1 && p->A771<0.0)
    msg = "A 771: the representative period must not be negative";
    else if(p->A770==1 && (p->A774<0 || p->A774>1))
    msg = "A 774: the long waves must be 0 (absorbing only) or 1 (bound long wave)";
    else if(p->A770==1 && (p->A748<0 || p->A748>1))
    msg = "A 748: the roller must be 0 (off) or 1 (on)";
    else if(p->A770==1 && !(p->A749>0.0))
    msg = "A 749: the roller slope must be positive";
    else if(p->A770==1 && !(p->A746>0.0))
    msg = "A 746: the Roelvink exponent must be positive";
    else if(p->A770==1 && p->A775!=1 && p->A775!=2)
    msg = "A 775: the advection of the wave groups must be 1 (first-order upwind) or 2 (second order)";
    else if(p->A770==1 && p->A747<0.0)
    msg = "A 747: the maximum H/h must not be negative";
    else if(p->A795<0.0 || p->A795>=1.0)
    msg = "A 795: the threshold of the spectral sparsity must be in [0,1)";
    else if(p->A795>0.0 && p->A700!=2)
    msg = "A 795: the spectral sparsity is for stationary runs (A 700 2)";
    else if(p->A796!=1 && p->A796!=2)
    msg = "A 796: the geographic advection must be 1 (first-order upwind) or 2 (second order)";
    else if(p->A796==2 && p->A770==1)
    msg = "A 796 2: not with the surfbeat model (its advection is set by A 775)";
    else if(p->A797!=0 && p->A797!=1)
    msg = "A 797: the sweeps with mesh refinement must be 0 (level by level) or 1 (composite)";
    else if(p->A715_k<1)
    msg = "A 715: the division of the sector directions must be at least 1";
    else if(p->A715_k>1 && p->A770==1)
    msg = "A 715: the fine direction sector is not available with the surfbeat model (A 770 1)";
    else if(p->A715_k>1 && p->A732==1 && p->A733==1)
    msg = "A 715: the DIA quadruplets (A 733 1) need uniform directions";
    else if(p->A738<1 || !(p->A739>=0.0))
    msg = "A 738, A 739: at least one source iteration, the tolerance must not be negative";
    else if(p->A798<1)
    msg = "A 798: at least one thread per rank";
    else if(p->A799!=0 && p->A799!=1)
    msg = "A 799: the convergence test must be 0 (change per iteration) or 1 (estimated distance to the solution)";

    if(msg!=nullptr)
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<msg<<endl<<endl;

        pgc->final(true);
    }
}

void seastate_f::environment(lexer *p, ghostcell *pgc)
{
    // bathymetry from the 2D grid; still water level F 60 (as SFLOW); a coupled host has set it
    if(!coupled)
    p->phimean = p->wd = p->F60;

    ILOOP
    JLOOP
    e->bed(i,j)=p->bed[IJ];

    pgc->gcsl_start4(p,e->bed,50);

    // bathymetry raster (A 790 1): the bed of every cell (halo included) from seastate-bathy.dat
    if(p->A790==1)
    {
    bathy = new seastate_bathy;
    std::string err;

        if(!bathy->read("seastate-bathy.dat",err))
        {
            if(p->mpirank==0)
            cout<<endl<<"SEASTATE input error  --  A 790 1: "<<err<<endl<<endl;

            pgc->final(true);
        }

        IMALOOP
        JMALOOP
        e->bed(i,j) = bathy->cell(p->XN[IP],p->XN[IP1],p->YN[JP],p->YN[JP1],p->wd-p->A705);

        if(p->mpirank==0)
        cout<<"SEASTATE bathymetry (A 790 1): seastate-bathy.dat, "<<bathy->nx<<" x "<<bathy->ny<<" nodes, spacing "<<bathy->dx<<" x "<<bathy->dy<<" m"<<endl;
    }

    IMALOOP
    JMALOOP
    {
    e->eta(i,j)=0.0;
    e->U(i,j)=0.0;
    e->V(i,j)=0.0;

    // prescribed current for stand-alone runs: U linear in x between xs and xe
    if(p->A720==1)
    {
    const double xs = p->A721_xs, xe = p->A721_xe;
    const double w = (xe>xs) ? std::min(std::max((p->XP[IP]-xs)/(xe-xs),0.0),1.0) : (p->XP[IP]>=xs ? 1.0 : 0.0);
    e->U(i,j) = (1.0-w)*p->A721_us + w*p->A721_ue;
    }
    const bool inside = i+p->origin_i>=0 && i+p->origin_i<p->gknox && j+p->origin_j>=0 && j+p->origin_j<p->gknoy;

    e->depth(i,j)=(inside && p->flagslice4[IJ]>0) ? std::max(p->wd - e->bed(i,j),0.0) : 0.0;
    e->wet(i,j)=(inside && p->flagslice4[IJ]>0 && e->depth(i,j)>=p->A705) ? 1 : 0;
    e->wet0(i,j)=e->wet(i,j);
    }
}

void seastate_f::storage(lexer *p, ghostcell *pgc)
{
    // surfbeat: one representative frequency, from the boundary spectrum on the grid A 701-703
    if(p->A770==1)
    {
    surfbeat_input(p,pgc);
    e->grid = new seastate_grid(1.0/trep,p->A703);
    }
    else
    {
    e->grid = new seastate_grid(p->A701,p->A702_fmin,p->A702_fmax,p->A703);

        // fine direction sector (A 715)
        if(p->A715_k>1 && e->grid->valid())
        e->grid->sector(p->A715_th1*3.14159265358979323846/180.0,p->A715_th2*3.14159265358979323846/180.0,p->A715_k);
    }

    // spectral sparsity (A 795): the ranges of directions are stored in one byte each
    if(p->A795>0.0 && e->grid->valid())
    {
    int nq[4] = {0,0,0,0};
    for(int m=0; m<e->grid->ndir; ++m)
    ++nq[e->grid->quad[m]];

        if(std::max(std::max(nq[0],nq[1]),std::max(nq[2],nq[3]))>250)
        {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  A 795: the spectral sparsity allows at most 250 directions per quadrant (A 703, A 715)"<<endl<<endl;

        pgc->final(true);
        }
    }

    if(!e->grid->valid())
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<e->grid->message()<<endl<<endl;

        pgc->final(true);
    }

    // block-sparse storage over the rank's index range including the ghost cells
    e->N = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,e->grid->nbin,p->A704);
    e->N->build(e->wet.V);

    // wave number and group velocity per frequency, same tiles
    e->kw = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,e->grid->nsig,p->A704);
    e->kw->build(e->wet.V);
    e->cg = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,e->grid->nsig,p->A704);
    e->cg->build(e->wet.V);

    // memory report
    int active=0;
    SLICELOOP4
    if(e->wet(i,j)==1)
    ++active;

    cells_active = pgc->globalsum(double(active));

    const double cells_alloc = pgc->globalsum(double(e->N->cells_allocated()));
    const double tiles_alloc = pgc->globalsum(double(e->N->tiles_allocated()));
    const double tiles_total = pgc->globalsum(double(e->N->tiles_total()));
    const double mb       = pgc->globalsum(double(e->N->bytes()))/1048576.0;
    const double mb_dense = pgc->globalsum(double(e->N->bytes_dense()))/1048576.0;
    const double mb_rank  = pgc->globalmax(double(e->N->bytes())/1048576.0);
    const double mb_kin   = pgc->globalsum(double(e->kw->bytes()+e->cg->bytes()))/1048576.0;

    if(p->mpirank==0)
    {
    const seastate_grid &g = *e->grid;

    cout<<endl<<"SEASTATE grid: "<<g.nsig<<" frequencies "<<fixed<<setprecision(3)<<g.fmin<<" - "<<g.fmax<<" Hz (ratio "<<setprecision(4)<<g.ratio<<"), "
        <<g.ndir<<" directions ("<<setprecision(2)<<g.dtheta*180.0/3.14159265358979323846<<" deg";
    if(!g.uniform)
    cout<<", fine sector "<<p->A715_th1<<" to "<<p->A715_th2<<" deg divided by "<<p->A715_k;
    cout<<"), "<<g.nbin<<" bins"<<endl;

    cout<<"SEASTATE storage: active cells "<<long(cells_active)<<" of "<<p->cellnumtot2D
        <<", allocated cells (incl. ghost cells) "<<long(cells_alloc)
        <<", tiles "<<long(tiles_alloc)<<" of "<<long(tiles_total)<<" ("<<p->A704<<" x "<<p->A704<<")"<<endl;

    cout<<"SEASTATE memory: N "<<setprecision(1)<<mb<<" MB float32 (dense: "<<mb_dense<<" MB), max per rank "<<mb_rank<<" MB; k and cg "<<mb_kin<<" MB"<<endl;
    if(p->A770==1)
    cout<<"SEASTATE mode: surfbeat, nonstationary, time step "<<(coupled ? "of the host" : "A 706")<<", refraction "<<p->A713<<", no frequency shift"<<endl;
    else
    cout<<"SEASTATE mode: "<<(p->A700==2 ? "stationary" : "nonstationary")<<", time step "<<p->A706<<" s, refraction "<<p->A713<<", frequency shift "<<p->A714<<endl;
    if(p->A795>0.0 || p->A796==2 || (p->A797==1 && p->G1>0))
    {
    cout<<"SEASTATE solver:";
    if(p->A795>0.0) cout<<" spectral sparsity (A 795, threshold "<<scientific<<setprecision(1)<<p->A795<<defaultfloat<<setprecision(6)<<" of the cell energy per bin)";
    if(p->A796==2) cout<<(p->A795>0.0 ? "," : "")<<" second-order geographic advection (A 796 2)";
    if(p->A797==1 && p->G1>0) cout<<(p->A795>0.0 || p->A796==2 ? "," : "")<<" composite sweep across the refinement levels (A 797 1)";
    cout<<endl;
    }
    if(p->A798>1 || (p->A700==2 && (p->A738>1 || p->A799==1)))
    {
    cout<<"SEASTATE solver (Phase 7b):";
    if(p->A798>1) cout<<" "<<p->A798<<" threads per rank (A 798)"<<(p->G1>0 || p->A770==1 ? ", not used with mesh refinement or surfbeat" : "");
    if(p->A700==2 && p->A738>1) cout<<(p->A798>1 ? "," : "")<<" up to "<<p->A738<<" source iterations per cell, tolerance "<<scientific<<setprecision(1)<<p->A739<<defaultfloat<<setprecision(6)<<" (A 738, A 739)";
    if(p->A700==2 && p->A799==1) cout<<(p->A798>1 || p->A738>1 ? "," : "")<<" convergence test on the estimated distance to the solution (A 799 1)";
    cout<<endl;
    }
    if(p->A720==1)
    cout<<"SEASTATE current: U "<<p->A721_us<<" m/s at x "<<p->A721_xs<<" m to "<<p->A721_ue<<" m/s at x "<<p->A721_xe<<" m"<<endl;
    cout<<endl;
    cout.unsetf(ios::floatfield);
    cout<<setprecision(6);
    }
}

void seastate_f::initial(lexer *p, ghostcell *pgc)
{
    if(p->A710==0)
    e->N->fill(0.0f);

    if(p->A710==1)
    initial_parametric(p,pgc);
}

void seastate_f::log_ini(lexer *p)
{
    if(p->mpirank!=0)
    return;

    mkdir("./REEF3D_SEASTATE_Log",0777);

    const char *path = "./REEF3D_SEASTATE_Log/REEF3D_SEASTATE_integral.dat";
    integral.open(path);

    integral<<"REEF3D::SEASTATE integral wave parameters"<<endl;
    integral<<"active cells: "<<long(cells_active)<<endl;
    integral<<"#iteration \t #simtime \t #E_tot [m^4] \t #Hs_max [m] \t #Hs_mean [m] \t #N_min \t #solver_iterations"<<endl;

    if(p->plog)
    p->plog->table_file(p,"seastate_integral","integral",path);
}

void seastate_f::log_step(lexer *p)
{
    if(p->mpirank!=0)
    return;

    integral<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<etot<<" \t "<<hsmax<<" \t "<<hsmean<<" \t "<<nmin<<" \t "<<iter_done<<endl;
}

void seastate_f::sources(lexer *p, ghostcell *pgc)
{
    seastate_source_param sp;

    sp.wind = (p->A730==1 || p->A730==2);
    sp.U10  = p->A731_u10;
    sp.wdir = p->A731_dir*3.14159265358979323846/180.0;
    sp.Alin = p->A734;

    sp.komen = (p->A732==1);
    sp.dia   = (p->A732==1 && p->A733==1);
    sp.limiter = p->A735;

    sp.breaking = (p->A740==1 || p->A740==2);
    sp.breaking_model = (p->A740==2) ? 2 : 1;
    sp.alpha    = p->A741_alpha;
    sp.gamma    = p->A741_gamma;
    sp.nroel    = p->A746;

    sp.friction = (p->A742==1);
    sp.Cb       = p->A743;

    sp.triads  = (p->A744==1);
    sp.alphaEB = p->A745;
    sp.cutfr   = p->A736_cutfr;
    sp.urcrit  = p->A736_urcrit;
    sp.urslim  = p->A736_urslim;
    sp.emax    = (p->A740==1 && p->A737==1);

    if(!sp.any())
    return;

    psrc = new seastate_source(*e->grid,sp);
    psolv->sources(psrc);

    if(p->mpirank==0)
    {
    cout<<"SEASTATE source terms:";
    if(sp.wind && p->A730==1)
    cout<<" wind U10 "<<sp.U10<<" m/s to "<<p->A731_dir<<" deg (Komen, linear growth "<<sp.Alin<<"),";
    if(sp.wind && p->A730==2)
    cout<<" wind field seastate-wind.dat (Komen, linear growth "<<sp.Alin<<"),";
    if(sp.komen)
    cout<<" whitecapping (Komen),";
    if(sp.dia)
    cout<<" quadruplets (DIA),";
    if(sp.komen)
    cout<<" action density limiter "<<sp.limiter<<",";
    if(sp.breaking && sp.breaking_model==1)
    cout<<" breaking (Battjes-Janssen, alpha "<<sp.alpha<<", gamma "<<sp.gamma<<"),";
    if(sp.breaking && sp.breaking_model==2)
    cout<<" breaking (Roelvink, alpha "<<sp.alpha<<", gamma "<<sp.gamma<<", n "<<sp.nroel<<"),";
    if(sp.friction)
    cout<<" bottom friction (JONSWAP, "<<sp.Cb<<" m^2/s^3),";
    if(sp.emax)
    cout<<" maximum energy (gamma d)^2/4,";
    if(sp.triads)
    cout<<" triads (LTA, alpha "<<sp.alphaEB<<", cutfr "<<sp.cutfr<<", urcrit "<<sp.urcrit<<", urslim "<<sp.urslim<<"),";
    cout<<endl;

    if(p->A700==2 && sp.komen)
    cout<<"SEASTATE stationary with the deep-water physics: pseudo time step "<<p->A706<<" s"<<endl;

    cout<<endl;
    }
}
