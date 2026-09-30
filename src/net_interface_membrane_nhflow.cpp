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

#include"net_interface.h"
#include"net_membrane.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"nhflow_membrane_beta.h"
#include<mpi.h>
#include<fstream>
#include<sstream>
#include<string>

// ctrl.txt: X 330 1 (and A 520 1). membrane.dat (read by rank 0, broadcast):
//
//   # comment
//   membrane box       x0 x1 y0 y1 z_bottom z_top     starts a new membrane (in 2D y0, y1 are ignored)
//   membrane cylinder  xc yc R z_bottom z_top
//   name        text                   optional label
//   resistance  R_n [R_t]              hydraulic resistance [m/s], leakage u_n = (dp/rho)/R_n; default 1e4, 0
//   thickness   delta                  half width of the smeared layer [m]; default 1.5 max(dx,dy,dz)
//   mesh        h                      target triangle edge length [m]; default min(dx,dy)
//   fill        dh                     initial inner water level above the outside level [m]; default 0
//   print       dt                     vtp output interval [s]; default none
//   poisson     0|1                    membrane mobility in the pressure Poisson equation; default 1
//                                      (0 only to demonstrate the splitting leakage of the projection)
//
// Parameter lines apply to the most recent 'membrane' line.

void net_interface::membrane_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    string content;
    int len=0;

    if(p->mpirank==0)
    {
        ifstream f("membrane.dat");

        if(!f)
        len=-1;
        else
        {
            stringstream ss;
            ss<<f.rdbuf();
            content = ss.str();
            len = (int)content.size();
        }
    }

    MPI_Bcast(&len,1,MPI_INT,0,pgc->mpi_comm);

    if(len<0)
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330: membrane.dat not found in the case directory !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    content.resize(len);
    if(len>0)
    MPI_Bcast(&content[0],len,MPI_CHAR,0,pgc->mpi_comm);

    // parse
    vector<membrane_param> mp;
    istringstream is(content);
    string line;
    int lineno=0;
    bool error=false;

    while(getline(is,line))
    {
        ++lineno;

        const size_t hash = line.find('#');
        if(hash!=string::npos)
        line = line.substr(0,hash);

        istringstream ls(line);
        string key;

        if(!(ls>>key))
        continue;

        if(key=="membrane")
        {
            string shape;
            membrane_param m;
            ls>>shape;

            if(shape=="box")
            {
                m.shape=1;
                if(!(ls>>m.x0>>m.x1>>m.y0>>m.y1>>m.zb>>m.zt))
                error=true;
            }
            else if(shape=="cylinder")
            {
                m.shape=2;
                if(!(ls>>m.xc>>m.yc>>m.R>>m.zb>>m.zt))
                error=true;
            }
            else
            error=true;

            m.name = "membrane"+to_string(mp.size());
            mp.push_back(m);
        }
        else if(mp.empty())
        error=true;
        else if(key=="name")
        {
            if(!(ls>>mp.back().name))
            error=true;
        }
        else if(key=="resistance")
        {
            if(!(ls>>mp.back().Rn))
            error=true;
            ls>>mp.back().Rt;
        }
        else if(key=="thickness")
        {
            if(!(ls>>mp.back().delta))
            error=true;
        }
        else if(key=="mesh")
        {
            if(!(ls>>mp.back().h))
            error=true;
        }
        else if(key=="fill")
        {
            if(!(ls>>mp.back().fill))
            error=true;
        }
        else if(key=="poisson")
        {
            if(!(ls>>mp.back().poisson))
            error=true;
        }
        else if(key=="print")
        {
            if(!(ls>>mp.back().printdt))
            error=true;
        }
        else
        error=true;

        if(error)
        {
            if(p->mpirank==0)
            cout<<"\n!!! membrane.dat, line "<<lineno<<": cannot read '"<<line<<"' !!!\n"<<endl;
            MPI_Abort(pgc->mpi_comm,1);
        }
    }

    if(mp.empty())
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330: no membrane defined in membrane.dat !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    // mobility field for the pressure Poisson equation, 1 away from membranes
    if(d->MBETA==nullptr)
    p->Darray(d->MBETA,p->imax*p->jmax*(p->kmax+2));

    if(d->MCHI==nullptr)
    p->Darray(d->MCHI,p->imax*p->jmax*(p->kmax+2));

    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    {
    d->MBETA[qn]=1.0;
    d->MCHI[qn]=0.0;
    }

    for(size_t m=0; m<mp.size(); ++m)
    {
        pmem.push_back(new net_membrane(m,mp[m]));
        pmem.back()->initialize_nhflow(p,d,pgc);
        pmem.back()->fill_nhflow(p,d,pgc);
    }
}

void net_interface::membrane_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha,
                                            double *UH, double *VH, double *WH, slice &WL)
{
    // 1. mobility field beta (also builds the membrane cell maps for this stage)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MBETA[qn]=1.0;

    for(auto m : pmem)
    m->mobility_nhflow(p,d,pgc,alpha);

    pgc->start4V(p,d->MBETA,1);

    // 2. static overpressure of the bag below its floor (prescribed pressure, see net_membrane)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MCHI[qn]=0.0;

    for(auto m : pmem)
    m->static_pressure_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);

    pgc->start4V(p,d->MCHI,1);

    // 3. implicit porous-jump forcing
    for(auto m : pmem)
    m->forcing_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
}

void net_interface::membrane_reaction_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, slice &WL, bool finalize)
{
    for(auto m : pmem)
    m->reaction_nhflow(p,d,pgc,alpha,WL,finalize);
}
