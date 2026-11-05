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

#include"nhflow_particle_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<cmath>
#include<limits>
#include<iostream>
#include<iomanip>

nhflow_particle_f::nhflow_particle_f(lexer *p, ghostcell *pgc) : gauss(0.0,1.0)
{
    static_assert(sizeof(nhflow_particle_data)==NF*sizeof(double),"nhflow_particle_data must be plain doubles");

    // local subdomain extent, model coordinates
    xloc_s = p->XN[0+marge];
    xloc_e = p->XN[p->knox+marge];

    if(p->j_dir==1)
    {
    yloc_s = p->YN[0+marge];
    yloc_e = p->YN[p->knoy+marge];
    }
    else
    {
    yloc_s = -std::numeric_limits<double>::max();
    yloc_e =  std::numeric_limits<double>::max();
    }

    // random numbers for the diffusion random walk, one stream per rank
    rng.seed(uint64_t(p->L11) + 104729ULL*uint64_t(p->mpirank));

    // wind for windage: L 41 overrides the NHFLOW wind forcing input A 571
    U10=0.0;
    double wdir=0.0;
    if(p->L41==1)
    {
    U10 = p->L41_u;
    wdir = p->L41_dir;
    }
    else if(p->A570>0)
    {
    U10 = p->A571_u;
    wdir = p->A571_dir;
    }
    cosw = cos(wdir*(PI/180.0));
    sinw = sin(wdir*(PI/180.0));

    setup_releases(p);

    printtime=0.0;
    printcount=0;
    numout=numout_global=0;

    if(p->mpirank==0)
    {
    mkdir("./REEF3D_NHFLOW_Particles",0777);

    logout.open("./REEF3D_NHFLOW_Particles/REEF3D-NHFLOW-Particles-Log.dat");
    logout<<"t \t released \t in_domain \t water_column \t surface \t stranded \t bed \t left_domain"<<endl;

    if(p->L62==2 || p->L62==3)
    {
    csvout.open("./REEF3D_NHFLOW_Particles/REEF3D-NHFLOW-Particles-Tracks.csv");
    csvout<<"t,id,release,x,y,z,u,v,w,state"<<endl;
    }

    long long total=0;
    for(auto &r : R)
    total+=r.num;

    cout<<"NHFLOW particles: "<<R.size()<<" releases, "<<total<<" particles, U10 windage: "<<U10<<" m/s"<<endl;
    }
}

nhflow_particle_f::~nhflow_particle_f()
{
    if(logout.is_open())
    logout.close();

    if(csvout.is_open())
    csvout.close();
}

void nhflow_particle_f::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // releases at t=0 appear in the initial output
    seed(p,d,pgc,p->simtime,true);

    for(auto &a : P)
    wetdry(p,d,a);

    print(p,d,pgc);
}

bool nhflow_particle_f::owned(lexer *p, double xp, double yp) const
{
    bool xin = (xp>=xloc_s && xp<xloc_e) || (p->nb4==-2 && xp>=xloc_s && xp<=xloc_e);

    if(p->j_dir==0)
    return xin;

    bool yin = (yp>=yloc_s && yp<yloc_e) || (p->nb2==-2 && yp>=yloc_s && yp<=yloc_e);

    return xin && yin;
}
