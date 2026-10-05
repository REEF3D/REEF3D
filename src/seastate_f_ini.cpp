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
    check_keys(p,pgc);

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

    environment(p,pgc);
    storage(p,pgc);
    kinematics(p,pgc,0.0);
    initial(p,pgc);
    boundary(p,pgc);

    // transport
    // nonstationary: 2 iterations on several ranks, so that the lagged halo values of the first
    // sweep are corrected (energy conservation across rank borders)
    iter_max = p->A707>0 ? p->A707 : (p->A700==2 ? 50 : (p->M10>1 ? 2 : 1));

    pex   = new seastate_exchange(p,e->grid->nbin,1);
    psolv = new seastate_implicit(p,e);

    sources(p,pgc);

    if(p->A700==1 && iter_max>1)
    {
    N0 = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,e->grid->nbin,p->A704);
    N0->build(e->wet.V);
    }

    pex->start(p,pgc,*e->N);

    parameters(p,pgc);

    pprint = new seastate_vtp(p,e,pgc);

    log_ini(p);

    // initial state
    pprint->start(p,e,pgc);
    log_step(p);
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
    else if(p->A711<0 || p->A711>2)
    msg = "A 711: boundary spectrum must be 0 (none), 1 (parametric) or 2 (SWAN spectrum file)";
    else if(p->A712_xm<0 || p->A712_xm>2 || p->A712_xp<0 || p->A712_xp>2 || p->A712_ym<0 || p->A712_ym>2 || p->A712_yp<0 || p->A712_yp>2)
    msg = "A 712: each side must be 0 (open), 1 (boundary spectrum) or 2 (zero gradient)";
    else if(p->A720!=0 && p->A720!=1)
    msg = "A 720: prescribed current must be 0 (none) or 1 (linear in x, A 721)";
    else if(p->A730!=0 && p->A730!=1)
    msg = "A 730: wind must be 0 (none) or 1 (uniform, A 731)";
    else if(p->A730==1 && !(p->A731_u10>=0.0))
    msg = "A 731: the wind speed must not be negative";
    else if(p->A730==1 && p->A732!=1)
    msg = "A 730 1: wind input needs the deep-water physics A 732 1 (Komen)";
    else if(p->A732!=0 && p->A732!=1)
    msg = "A 732: deep-water physics must be 0 (off) or 1 (Komen)";
    else if(p->A733!=0 && p->A733!=1)
    msg = "A 733: quadruplets must be 0 (off) or 1 (DIA)";
    else if(!(p->A734>=0.0))
    msg = "A 734: the linear growth coefficient must not be negative";
    else if(!(p->A735>=0.0))
    msg = "A 735: the limiter coefficient must not be negative";
    else if(p->A740!=0 && p->A740!=1)
    msg = "A 740: depth-induced breaking must be 0 (off) or 1 (Battjes-Janssen)";
    else if(p->A740==1 && !(p->A741_alpha>0.0 && p->A741_gamma>0.0))
    msg = "A 741: the breaking coefficients alpha and gamma must be positive";
    else if(p->A742!=0 && p->A742!=1)
    msg = "A 742: bottom friction must be 0 (off) or 1 (JONSWAP)";
    else if(!(p->A743>=0.0))
    msg = "A 743: the friction coefficient must not be negative";
    else if(p->A744!=0 && p->A744!=1)
    msg = "A 744: triads must be 0 (off) or 1 (LTA)";
    else if(!(p->A745>=0.0))
    msg = "A 745: the triad coefficient must not be negative";

    if(msg!=nullptr)
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<msg<<endl<<endl;

        pgc->final(true);
    }
}

void seastate_f::environment(lexer *p, ghostcell *pgc)
{
    // bathymetry from the 2D grid; still water level F 60 (as SFLOW)
    p->phimean = p->wd = p->F60;

    ILOOP
    JLOOP
    e->bed(i,j)=p->bed[IJ];

    pgc->gcsl_start4(p,e->bed,50);

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
    }
}

void seastate_f::storage(lexer *p, ghostcell *pgc)
{
    e->grid = new seastate_grid(p->A701,p->A702_fmin,p->A702_fmax,p->A703);

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
        <<g.ndir<<" directions ("<<setprecision(2)<<g.dtheta*180.0/3.14159265358979323846<<" deg), "<<g.nbin<<" bins"<<endl;

    cout<<"SEASTATE storage: active cells "<<long(cells_active)<<" of "<<p->cellnumtot2D
        <<", allocated cells (incl. ghost cells) "<<long(cells_alloc)
        <<", tiles "<<long(tiles_alloc)<<" of "<<long(tiles_total)<<" ("<<p->A704<<" x "<<p->A704<<")"<<endl;

    cout<<"SEASTATE memory: N "<<setprecision(1)<<mb<<" MB float32 (dense: "<<mb_dense<<" MB), max per rank "<<mb_rank<<" MB; k and cg "<<mb_kin<<" MB"<<endl;
    cout<<"SEASTATE mode: "<<(p->A700==2 ? "stationary" : "nonstationary")<<", time step "<<p->A706<<" s, refraction "<<p->A713<<", frequency shift "<<p->A714<<endl;
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

    sp.wind = (p->A730==1);
    sp.U10  = p->A731_u10;
    sp.wdir = p->A731_dir*3.14159265358979323846/180.0;
    sp.Alin = p->A734;

    sp.komen = (p->A732==1);
    sp.dia   = (p->A732==1 && p->A733==1);
    sp.limiter = p->A735;

    sp.breaking = (p->A740==1);
    sp.alpha    = p->A741_alpha;
    sp.gamma    = p->A741_gamma;

    sp.friction = (p->A742==1);
    sp.Cb       = p->A743;

    sp.triads  = (p->A744==1);
    sp.alphaEB = p->A745;

    if(!sp.any())
    return;

    psrc = new seastate_source(*e->grid,sp);
    psolv->sources(psrc);

    if(p->mpirank==0)
    {
    cout<<"SEASTATE source terms:";
    if(sp.wind)
    cout<<" wind U10 "<<sp.U10<<" m/s to "<<p->A731_dir<<" deg (Komen, linear growth "<<sp.Alin<<"),";
    if(sp.komen)
    cout<<" whitecapping (Komen),";
    if(sp.dia)
    cout<<" quadruplets (DIA),";
    if(sp.komen)
    cout<<" action density limiter "<<sp.limiter<<",";
    if(sp.breaking)
    cout<<" breaking (Battjes-Janssen, alpha "<<sp.alpha<<", gamma "<<sp.gamma<<"),";
    if(sp.friction)
    cout<<" bottom friction (JONSWAP, "<<sp.Cb<<" m^2/s^3),";
    if(sp.triads)
    cout<<" triads (LTA, alpha "<<sp.alphaEB<<"),";
    cout<<endl;

    if(p->A700==2 && sp.komen)
    cout<<"SEASTATE stationary with the deep-water physics: pseudo time step "<<p->A706<<" s"<<endl;

    cout<<endl;
    }
}
