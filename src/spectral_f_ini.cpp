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

#include"spectral_f.h"
#include"fdm_spectral.h"
#include"spectral_grid.h"
#include"spectral_store.h"
#include"spectral_vtp.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include<algorithm>
#include<iostream>
#include<iomanip>
#include<sys/stat.h>
#include<sys/types.h>

void spectral_f::ini(lexer *p, ghostcell *pgc)
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
    initial(p,pgc);
    parameters(p,pgc);

    pprint = new spectral_vtp(p,e,pgc);

    log_ini(p);

    // initial state
    pprint->start(p,e,pgc);
    log_step(p);
}

void spectral_f::check_keys(lexer *p, ghostcell *pgc)
{
    const char *msg = nullptr;

    if(p->A704<1)
    msg = "A 704: the tile size must be at least 1";
    else if(p->A705<0.0)
    msg = "A 705: the minimum water depth must not be negative";
    else if(!(p->A706>0.0))
    msg = "A 706: the time step must be positive";
    else if(p->A710!=0 && p->A710!=1)
    msg = "A 710: initial spectrum must be 0 (zero) or 1 (parametric)";

    if(msg!=nullptr)
    {
        if(p->mpirank==0)
        cout<<endl<<"Spectral input error  --  "<<msg<<endl<<endl;

        pgc->final(true);
    }
}

void spectral_f::environment(lexer *p, ghostcell *pgc)
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
    e->depth(i,j)=p->flagslice4[IJ]>0 ? std::max(p->wd - e->bed(i,j),0.0) : 0.0;
    e->wet(i,j)=(p->flagslice4[IJ]>0 && e->depth(i,j)>=p->A705) ? 1 : 0;
    }
}

void spectral_f::storage(lexer *p, ghostcell *pgc)
{
    e->grid = new spectral_grid(p->A701,p->A702_fmin,p->A702_fmax,p->A703);

    if(!e->grid->valid())
    {
        if(p->mpirank==0)
        cout<<endl<<"Spectral input error  --  "<<e->grid->message()<<endl<<endl;

        pgc->final(true);
    }

    // block-sparse storage over the rank's index range including the ghost cells
    e->N = new spectral_store(p->imin,p->jmin,p->imax,p->jmax,e->grid->nbin,p->A704);
    e->N->build(e->wet.V);

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

    if(p->mpirank==0)
    {
    const spectral_grid &g = *e->grid;

    cout<<endl<<"Spectral grid: "<<g.nsig<<" frequencies "<<fixed<<setprecision(3)<<g.fmin<<" - "<<g.fmax<<" Hz (ratio "<<setprecision(4)<<g.ratio<<"), "
        <<g.ndir<<" directions ("<<setprecision(2)<<g.dtheta*180.0/3.14159265358979323846<<" deg), "<<g.nbin<<" bins"<<endl;

    cout<<"Spectral storage: active cells "<<long(cells_active)<<" of "<<p->cellnumtot2D
        <<", allocated cells (incl. ghost cells) "<<long(cells_alloc)
        <<", tiles "<<long(tiles_alloc)<<" of "<<long(tiles_total)<<" ("<<p->A704<<" x "<<p->A704<<")"<<endl;

    cout<<"Spectral memory: N "<<setprecision(1)<<mb<<" MB float32 (dense: "<<mb_dense<<" MB), max per rank "<<mb_rank<<" MB"<<endl<<endl;
    cout.unsetf(ios::floatfield);
    cout<<setprecision(6);
    }
}

void spectral_f::initial(lexer *p, ghostcell *pgc)
{
    if(p->A710==0)
    e->N->fill(0.0f);

    if(p->A710==1)
    initial_parametric(p,pgc);
}

void spectral_f::log_ini(lexer *p)
{
    if(p->mpirank!=0)
    return;

    mkdir("./REEF3D_SPECTRAL_Log",0777);

    const char *path = "./REEF3D_SPECTRAL_Log/REEF3D_SPECTRAL_integral.dat";
    integral.open(path);

    integral<<"REEF3D::Spectral integral wave parameters"<<endl;
    integral<<"active cells: "<<long(cells_active)<<endl;
    integral<<"#iteration \t #simtime \t #E_tot [m^4] \t #Hs_max [m] \t #Hs_mean [m]"<<endl;

    if(p->plog)
    p->plog->table_file(p,"spectral_integral","integral",path);
}

void spectral_f::log_step(lexer *p)
{
    if(p->mpirank!=0)
    return;

    integral<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<etot<<" \t "<<hsmax<<" \t "<<hsmean<<endl;
}
