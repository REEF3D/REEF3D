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

#include"fem_coupling.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdio>
#include<algorithm>
#include<sys/stat.h>
#include<sys/types.h>

fem_coupling::fem_coupling(lexer *p, ghostcell *pgc) : surf_version(-1), dxmin(0.0), rho_w(1000.0), printtime(0.0), printcount(0), starttime(0.0)
{
    // rank 0 reads fem.dat, everybody parses the broadcast content
    std::string content;
    int len = 0;

    if(p->mpirank==0)
    {
        std::ifstream f("fem.dat");
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
        std::cout<<"\n!!! Z 30: fem.dat not found in the case directory !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    content.resize(len);
    if(len>0)
    MPI_Bcast(&content[0],len,MPI_CHAR,0,pgc->mpi_comm);

    try
    {
        fs.set_gravity(fem_solid::Vec3(p->W20,p->W21,p->W22));
        if(p->j_dir==0)
        fs.set_plane_strain(true);

        std::istringstream is(content);
        fs.read(is);
        fs.build();
    }
    catch(std::exception& e)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! "<<e.what()<<" !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    rho_w = p->W1;

    // smallest fluid cell size: spacing of the Lagrangian points
    double h = 1.0e20;
    for(int ii=0; ii<p->knox; ++ii)
    h = std::min(h,p->DXN[ii+marge]);
    if(p->j_dir==1)
    for(int jj=0; jj<p->knoy; ++jj)
    h = std::min(h,p->DYN[jj+marge]);
    for(int kk=0; kk<p->knoz; ++kk)
    h = std::min(h,p->DZN[kk+marge]);
    dxmin = pgc->globalmin(h);

    ini_points(p,pgc);

    outdir = "./REEF3D_FEM";

    if(p->mpirank==0)
    {
        mkdir(outdir.c_str(),0777);

        fs.info(std::cout);
        std::cout<<"FEM: "<<pts.size()<<" Lagrangian points on the surface, critical solid time step "<<fs.critical_dt()<<" s"<<std::endl;

        std::ofstream lg((outdir+"/REEF3D_FEM_log.dat").c_str());
        lg<<"# time  substeps  intact_elements  eroded  debris  F_fluid_x F_fluid_y F_fluid_z  F_support_x F_support_y F_support_z  max_vonMises  max_displacement  kinetic_energy  dissipated_energy  (SI units, support force = force of the structure on its supports)\n";

        for(const fem_solid::monitor& mo : fs.monitors())
        {
            std::ofstream ts((outdir+"/REEF3D_FEM_monitor_"+mo.name+".dat").c_str());
            const fem_solid::Vec3& X = fs.ref_pos(mo.node);
            ts<<"# node at "<<X(0)<<" "<<X(1)<<" "<<X(2)<<"\n# time  dx dy dz  vx vy vz\n";
        }
    }

    print(p);
}

fem_coupling::~fem_coupling()
{
}

void fem_coupling::ini_points(lexer *p, ghostcell *pgc)
{
    // Lagrangian points: every surface face is split into n1 x n2 points with
    // a spacing below the smallest fluid cell. In 2D (j_dir 0) the faces
    // normal to y are the out-of-plane boundaries and carry no points.
    (void)pgc;
    pts.clear();

    const std::vector<fem_solid::face>& F = fs.surface();
    const int npf = fs.coupling().points_per_face;

    for(int f=0; f<(int)F.size(); ++f)
    {
        const fem_solid::Vec3& x0 = fs.ref_pos(F[f].n[0]);
        const fem_solid::Vec3& x1 = fs.ref_pos(F[f].n[1]);
        const fem_solid::Vec3& x3 = fs.ref_pos(F[f].n[3]);

        const fem_solid::Vec3 nrm = (fs.ref_pos(F[f].n[2])-x0).cross(x3-x1);
        if(p->j_dir==0 && std::fabs(nrm(1))>0.5*nrm.norm())
        continue;

        int n1 = npf>0 ? npf : std::max(1,(int)std::ceil((x1-x0).norm()/(0.9*dxmin)));
        int n2 = npf>0 ? npf : std::max(1,(int)std::ceil((x3-x0).norm()/(0.9*dxmin)));

        // 2D: one point across the out-of-plane direction
        if(p->j_dir==0)
        {
            const fem_solid::Vec3 e1 = x1-x0, e2 = x3-x0;
            if(std::fabs(e1(1))>0.5*e1.norm()) n1 = 1;
            if(std::fabs(e2(1))>0.5*e2.norm()) n2 = 1;
        }

        for(int a=0; a<n1; ++a)
        for(int b=0; b<n2; ++b)
        {
            lpoint q;
            q.face = f;
            q.s = (a+0.5)/double(n1);
            q.t = (b+0.5)/double(n2);
            q.frac = 1.0/double(n1*n2);
            pts.push_back(q);
        }
    }

    surf_version = fs.surface_version();

    buf.assign(8*(pts.size()+fs.debris().size()) + 2*size_t(fs.nnode()),0.0);
    fdeb.assign(fs.debris().size(),fem_solid::Vec3::Zero());
}

void fem_coupling::point_state(int q, fem_solid::Vec3& xp, fem_solid::Vec3& vp, fem_solid::Vec3& n, double& A) const
{
    const lpoint& L = pts[q];
    const fem_solid::face& F = fs.surface()[L.face];

    // node order n0 -> n1 -> n2 -> n3 around the face, s along n0-n1, t along n0-n3
    const double w[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};

    xp.setZero();
    vp.setZero();
    for(int a=0; a<4; ++a)
    {
        xp += w[a]*fs.pos(F.n[a]);
        vp += w[a]*fs.vel(F.n[a]);
    }

    const fem_solid::Vec3 c = (fs.pos(F.n[2])-fs.pos(F.n[0])).cross(fs.pos(F.n[3])-fs.pos(F.n[1]));
    const double cn = c.norm();
    n = cn>0.0 ? fem_solid::Vec3(c/cn) : fem_solid::Vec3::Zero();
    A = 0.5*cn*L.frac;
}

double fem_coupling::kernel(double r) const
{
    // 3-point kernel of Roma et al. (1999), as in the FSI strips
    r = std::fabs(r);
    if(r<=0.5)
    return (1.0 + std::sqrt(1.0 - 3.0*r*r))/3.0;
    if(r<=1.5)
    return (5.0 - 3.0*r - std::sqrt(1.0 - 3.0*(1.0-r)*(1.0-r)))/6.0;
    return 0.0;
}

void fem_coupling::print(lexer *p)
{
    if(p->mpirank!=0)
    return;

    const double dtp = fs.coupling().print_dt>0.0 ? fs.coupling().print_dt : p->Z31;

    if(dtp>0.0 && fs.time()>=printtime-1.0e-12)
    {
        char name[256];
        std::snprintf(name,sizeof(name),"%s/REEF3D-FEM-%08d.vtu",outdir.c_str(),printcount);
        try
        {
            fs.write_vtu(name);
        }
        catch(std::exception& e)
        {
            std::cout<<e.what()<<std::endl;
        }
        ++printcount;
        printtime += dtp;
    }
}
