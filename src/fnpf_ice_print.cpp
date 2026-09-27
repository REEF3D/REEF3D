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

#include"fnpf_ice.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<cstdio>
#include<iomanip>

void fnpf_ice::print_ini(lexer *p)
{
    printtime = 0.0;
    printcount = 0;

    if(p->mpirank!=0)
    return;

    mkdir("./REEF3D_FNPF_ICE",0777);

    logout.open("./REEF3D_FNPF_ICE/REEF3D-FNPF-ICE-floes.dat");
    logout<<"# FNPF ice floes: time, id, centre of mass, Euler angles ZYX [deg], velocity, hydrodynamic force, contact force (planar), loaded area"<<endl;
    logout<<"# t id x y z roll pitch yaw u v w Fx_hyd Fy_hyd Fz_hyd Fx_c Fy_c A_wet sigma_f"<<endl;

    if(breakflag>0)
    {
    breakout.open("./REEF3D_FNPF_ICE/REEF3D-FNPF-ICE-breaking.dat");
    breakout<<"# ice breaking events: mech 1 flexural (value, limit: stress, strength [Pa]), 2 contact splitting (force, splitting load [N]),"<<endl;
    breakout<<"# 3 spalling (crushed depth, spall length [m]); new_id -2: chip cleared as rubble"<<endl;
    breakout<<"# t parent_id new_id mech value limit area_parent_piece area_new_piece"<<endl;
    }
    
    if(p->A398>0.0)
    {
    crushout.open("./REEF3D_FNPF_ICE/REEF3D-FNPF-ICE-crushing.dat");
    crushout<<"# ice crushing contacts, every time step with a capped or penetrating pair (A 398 sigma_c = "<<p->A398<<" Pa)"<<endl;
    crushout<<"# id_b < 0: domain wall. depth: crushed depth (overlap), volume = h*width*depth"<<endl;
    crushout<<"# t id_a id_b Fn Fcap depth width volume"<<endl;
    }

    if(nobst>0)
    {
    obstout.open("./REEF3D_FNPF_ICE/REEF3D-FNPF-ICE-obstacle-loads.dat");
    obstout<<"# ice loads on the obstacles (contact force, planar), every time step"<<endl;
    obstout<<"# t";
    for(auto &fl : floe)
    if(fl.type==1)
    obstout<<"  Fx_"<<fl.id<<" Fy_"<<fl.id;
    obstout<<endl;
    }
}

void fnpf_ice::print(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->mpirank!=0)
    return;

    const double t = (p->count==0) ? p->simtime : p->simtime + p->dt;

    // obstacle loads every step
    if(obstout.is_open() && p->count>0)
    {
    obstout<<setprecision(9)<<t;
    for(auto &fl : floe)
    if(fl.type==1)
    obstout<<" "<<fl.Fc[0]<<" "<<fl.Fc[1];
    obstout<<endl;
    }

    // print interval: A 387 in s, otherwise the fsf vtp settings (P 180, P 181, P 182)
    int doprint=0;

    if(p->count==0)
    doprint=1;

    else if(p->A387>0.0)
    doprint = (t>=printtime);

    else if(p->P182>0.0)
    doprint = (t>=printtime);

    else if(p->P181>0)
    doprint = (p->count%p->P181==0);

    else
    doprint = (p->count%100==0);

    if(doprint==0)
    return;

    const double interval = (p->A387>0.0) ? p->A387 : p->P182;
    if(interval>0.0)
    while(printtime<=t)
    printtime += interval;

    print_vtp(p,printcount);
    print_log(p);
    ++printcount;
}

void fnpf_ice::print_log(lexer *p)
{
    const double t = (p->count==0) ? p->simtime : p->simtime + p->dt;
    const double deg = 180.0/PI;

    for(auto &fl : floe)
    {
        if(fl.type!=0)
        continue;

        double roll,pitch,yaw;
        quat_to_euler(fl.q,roll,pitch,yaw);

        logout<<setprecision(9)<<t<<" "<<fl.id<<" "
              <<fl.x[0]<<" "<<fl.x[1]<<" "<<fl.x[2]<<" "
              <<roll*deg<<" "<<pitch*deg<<" "<<yaw*deg<<" "
              <<fl.v[0]<<" "<<fl.v[1]<<" "<<fl.v[2]<<" "
              <<fl.F[0]<<" "<<fl.F[1]<<" "<<fl.F[2]<<" "
              <<fl.Fc[0]<<" "<<fl.Fc[1]<<" "<<fl.Awet<<" "<<fl.sigf<<endl;
    }
}

void fnpf_ice::print_vtp(lexer *p, int num)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_FNPF_ICE/REEF3D-FNPF-ICE-%08i.vtp",num);

    ofstream out(name);

    // removed floes are not drawn, rubble (type 3) as floes
    vector<fnpf_ice_floe> vis;
    for(size_t n=0; n<floe.size(); ++n)
    if(floe[n].type!=2)
    vis.push_back(floe[n]);

    // prisms: bottom and top polygons plus side quads; obstacles drawn as columns through the surface
    double hob=0.0;
    for(auto &fl : vis)
    if(fl.type==0)
    hob = MAX(hob,fl.h);
    if(hob<=0.0)
    hob=1.0;

    int npts=0, npoly=0, nconn=0;
    for(auto &fl : vis)
    {
    const int nv = int(fl.bx.size());
    npts  += 2*nv;
    npoly += 2 + nv;
    nconn += 2*nv + 4*nv;
    }

    out<<"<?xml version=\"1.0\"?>"<<endl;
    out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    out<<"<PolyData>"<<endl;
    out<<"<Piece NumberOfPoints=\""<<npts<<"\" NumberOfPolys=\""<<npoly<<"\">"<<endl;

    out<<"<Points>"<<endl;
    out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    for(auto &fl : vis)
    {
        double R[3][3];
        quat_to_matrix(fl.q,R);
        const int nv = int(fl.bx.size());

        for(int side=0; side<2; ++side)
        for(int q=0; q<nv; ++q)
        {
            if(fl.type==0 || fl.type==3)
            {
            const double zb = side==0 ? -0.5*fl.h : 0.5*fl.h;
            out<<fl.x[0] + R[0][0]*fl.bx[q] + R[0][1]*fl.by[q] + R[0][2]*zb<<" "
               <<fl.x[1] + R[1][0]*fl.bx[q] + R[1][1]*fl.by[q] + R[1][2]*zb<<" "
               <<fl.x[2] + R[2][0]*fl.bx[q] + R[2][1]*fl.by[q] + R[2][2]*zb<<endl;
            }
            else
            out<<fl.x[0]+fl.bx[q]<<" "<<fl.x[1]+fl.by[q]<<" "<<wd + (side==0 ? -2.0*hob : 2.0*hob)<<endl;
        }
    }
    out<<"</DataArray>"<<endl;
    out<<"</Points>"<<endl;

    out<<"<CellData Scalars=\"id\">"<<endl;
    out<<"<DataArray type=\"Int32\" Name=\"id\" format=\"ascii\">"<<endl;
    for(auto &fl : vis)
    for(size_t q=0; q<fl.bx.size()+2; ++q)
    out<<fl.id<<endl;
    out<<"</DataArray>"<<endl;
    out<<"<DataArray type=\"Int32\" Name=\"type\" format=\"ascii\">"<<endl;
    for(auto &fl : vis)
    for(size_t q=0; q<fl.bx.size()+2; ++q)
    out<<fl.type<<endl;
    out<<"</DataArray>"<<endl;
    out<<"<DataArray type=\"Float32\" Name=\"speed\" format=\"ascii\">"<<endl;
    for(auto &fl : vis)
    for(size_t q=0; q<fl.bx.size()+2; ++q)
    out<<sqrt(fl.v[0]*fl.v[0]+fl.v[1]*fl.v[1]+fl.v[2]*fl.v[2])<<endl;
    out<<"</DataArray>"<<endl;
    out<<"<DataArray type=\"Float32\" Name=\"contact_force\" format=\"ascii\">"<<endl;
    for(auto &fl : vis)
    for(size_t q=0; q<fl.bx.size()+2; ++q)
    out<<sqrt(fl.Fc[0]*fl.Fc[0]+fl.Fc[1]*fl.Fc[1])<<endl;
    out<<"</DataArray>"<<endl;
    out<<"</CellData>"<<endl;

    out<<"<Polys>"<<endl;
    out<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<endl;
    int base=0;
    for(auto &fl : vis)
    {
        const int nv = int(fl.bx.size());
        // bottom (reversed, facing down), top
        for(int q=nv-1; q>=0; --q) out<<base+q<<" ";
        out<<endl;
        for(int q=0; q<nv; ++q) out<<base+nv+q<<" ";
        out<<endl;
        for(int q=0; q<nv; ++q)
        {
        const int q2=(q+1)%nv;
        out<<base+q<<" "<<base+q2<<" "<<base+nv+q2<<" "<<base+nv+q<<endl;
        }
        base += 2*nv;
    }
    out<<"</DataArray>"<<endl;
    out<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<endl;
    int off=0;
    for(auto &fl : vis)
    {
        const int nv = int(fl.bx.size());
        off+=nv; out<<off<<endl;
        off+=nv; out<<off<<endl;
        for(int q=0; q<nv; ++q)
        {
        off+=4;
        out<<off<<endl;
        }
    }
    out<<"</DataArray>"<<endl;
    out<<"</Polys>"<<endl;

    out<<"</Piece>"<<endl;
    out<<"</PolyData>"<<endl;
    out<<"</VTKFile>"<<endl;

    out.close();

    (void)nconn;
}
