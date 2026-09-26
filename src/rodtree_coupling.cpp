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

#include"rodtree_coupling.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdio>
#include<sys/stat.h>
#include<sys/types.h>

rodtree_coupling::rodtree_coupling(lexer *p, ghostcell *pgc) : npts(0), tprev(0.0), have_prev(false), printtime(0.0), printcount(0)
{
    // rank 0 reads rodtree.dat, everybody parses the broadcast content
    std::string content;
    int len = 0;

    if(p->mpirank==0)
    {
        std::ifstream f("rodtree.dat");
        if(!f)
        len = -1;
        else
        {
            std::stringstream ss;
            ss<<f.rdbuf();
            content = ss.str();
            len = (int)content.size();
        }
    }

    MPI_Bcast(&len,1,MPI_INT,0,pgc->mpi_comm);

    if(len<0)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! Z 20: rodtree.dat not found in the case directory !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    content.resize(len);
    if(len>0)
    MPI_Bcast(&content[0],len,MPI_CHAR,0,pgc->mpi_comm);

    try
    {
        std::istringstream is(content);
        rt.read(is);
        rt.finalize_setup();
    }
    catch(std::exception& e)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! "<<e.what()<<" !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    rho = p->W1;
    rt.set_fluid_density(rho);
    rt.set_gravity(Eigen::Vector3d(p->W20,p->W21,p->W22));

    ini_points(p,pgc);

    ufprev.assign(rt.nelem(),Eigen::Vector3d::Zero());

    outdir = "./REEF3D_RodTree";

    if(p->mpirank==0)
    {
        mkdir(outdir.c_str(),0777);

        std::cout<<"RodTree: "<<rt.ncolony()<<" colonies, "<<rt.nelem()<<" elements, "<<npts<<" Lagrangian points, "
                 <<(rt.get_integrator()==0 ? "implicit" : "explicit")<<" integrator";
        if(rt.get_integrator()==0)
        std::cout<<" ("<<rt.get_substeps()<<" substeps)";
        std::cout<<", reaction "<<(rt.get_reaction_mode()==0 ? "full" : rt.get_reaction_mode()==1 ? "drag" : "none")<<std::endl;

        for(int c=0; c<rt.ncolony(); ++c)
        {
            std::ofstream ts((outdir + "/REEF3D_RodTree_" + rt.colony(c).name + ".dat").c_str());
            ts<<"# time  tip_dx tip_dy tip_dz  base_Fx base_Fy base_Fz  hydro_Fx hydro_Fy hydro_Fz  submerged_elements  polyp_extension  (SI units; base force = force of the colony on the bed)\n";
        }
    }
}

rodtree_coupling::~rodtree_coupling()
{
}

void rodtree_coupling::ini_points(lexer *p, ghostcell *pgc)
{
    // Lagrangian spacing <= smallest horizontal cell size
    double h = 1.0e20;

    for(int ii=0; ii<p->knox; ++ii)
    h = std::min(h,p->DXN[ii+marge]);

    if(p->j_dir==1)
    for(int jj=0; jj<p->knoy; ++jj)
    h = std::min(h,p->DYN[jj+marge]);

    h = pgc->globalmin(h);

    const int ne = rt.nelem();
    nq.resize(ne);
    first.resize(ne);
    npts = 0;

    for(int e=0; e<ne; ++e)
    {
        nq[e] = std::max(1,std::min(64,(int)std::ceil(rt.elem(e).l/h - 1.0e-9)));
        first[e] = npts;
        npts += nq[e];
    }

    buf.assign(4*npts,0.0);
}

void rodtree_coupling::reduce_samples(ghostcell *pgc)
{
    if(npts>0)
    MPI_Allreduce(MPI_IN_PLACE,buf.data(),4*npts,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
}

void rodtree_coupling::apply_samples()
{
    for(int e=0; e<rt.nelem(); ++e)
    {
        Eigen::Vector3d us = Eigen::Vector3d::Zero();
        int nin = 0;

        for(int q=0; q<nq[e]; ++q)
        {
            const double *b = &buf[4*(first[e]+q)];
            if(b[3]>0.5)
            {
                us += Eigen::Vector3d(b[0],b[1],b[2])/b[3];
                ++nin;
            }
        }

        rodtree::element& E = rt.elem(e);
        E.chi = double(nin)/double(nq[e]);
        E.uf = nin>0 ? Eigen::Vector3d(us/double(nin)) : Eigen::Vector3d::Zero();
    }
}

void rodtree_coupling::finish_step(lexer *p, ghostcell *pgc, bool finalize)
{
    if(!finalize)
    return;

    // local fluid acceleration at the elements (backward difference over the step)
    for(int e=0; e<rt.nelem(); ++e)
    {
        rodtree::element& E = rt.elem(e);
        E.af = (have_prev && p->dt>0.0) ? Eigen::Vector3d((E.uf - ufprev[e])/p->dt) : Eigen::Vector3d::Zero();
        ufprev[e] = E.uf;
    }
    have_prev = true;

    try
    {
        rt.advance(p->dt);
    }
    catch(std::exception& e)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! "<<e.what()<<" !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    print(p);
}

void rodtree_coupling::print(lexer *p)
{
    if(p->mpirank!=0)
    return;

    for(int c=0; c<rt.ncolony(); ++c)
    {
        Eigen::Vector3d dx = rt.tip_displacement(c), fb = rt.base_force(c), fh = rt.hydro_force(c);
        int nsub = 0;
        for(int e : rt.colony(c).elements)
        if(rt.elem(e).chi>0.0) ++nsub;

        std::ofstream ts((outdir + "/REEF3D_RodTree_" + rt.colony(c).name + ".dat").c_str(),std::ios::app);
        ts<<std::setprecision(9)<<rt.time()<<" "<<dx(0)<<" "<<dx(1)<<" "<<dx(2)<<" "
          <<fb(0)<<" "<<fb(1)<<" "<<fb(2)<<" "<<fh(0)<<" "<<fh(1)<<" "<<fh(2)<<" "<<nsub<<" "<<rt.polyp_extension(c)<<"\n";
    }

    if(p->Z21>0.0 && rt.time()>=printtime-1.0e-12)
    {
        char name[256];
        std::snprintf(name,sizeof(name),"%s/REEF3D-RodTree-%08d.vtp",outdir.c_str(),printcount);
        rt.write_vtp(name);
        ++printcount;
        printtime += p->Z21;
    }
}

double rodtree_coupling::kernel(double r) const
{
    r = std::fabs(r);
    return r<2.0 ? 0.25*(1.0 + std::cos(0.5*PI*r)) : 0.0;
}

double rodtree_coupling::self_weight(double frac) const
{
    // sum of squared 1D kernel weights on a uniform grid for a point at
    // offset frac (cell units) from a cell centre: fluid response of one
    // point to its own spread force is F*S/(rho*Vcell)
    double s = 0.0;
    for(int kk=-3; kk<=3; ++kk)
    {
        double w = kernel(double(kk) - frac);
        s += w*w;
    }
    return s;
}
