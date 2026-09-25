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
#include<mpi.h>
#include<cmath>
#include<cstdio>
#include<cstdint>
#include<iomanip>

// output every L 61 seconds of simulation time
// L 62 = 1: VTP (ParaView), 2: CSV tracks, 3: both
// all particles are gathered on rank 0, positions and velocities in world coordinates

void nhflow_particle_f::print(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double eps = 1.0e-10;

    if(p->simtime < printtime - eps)
    return;

    std::vector<nhflow_particle_data> all;
    gather(p,pgc,all);

    if(p->mpirank==0)
    {
        if(p->L62==1 || p->L62==3)
        print_vtp(p,all);

        if(p->L62==2 || p->L62==3)
        print_csv(p,all);
    }

    print_log(p,pgc);

    ++printcount;

    if(p->L61>0.0)
    while(printtime <= p->simtime + eps)
    printtime += p->L61;
}

void nhflow_particle_f::gather(lexer *p, ghostcell *pgc, std::vector<nhflow_particle_data> &all)
{
    int num = int(P.size());
    std::vector<int> nums(p->mpi_size,0), counts(p->mpi_size,0), displ(p->mpi_size,0);

    MPI_Gather(&num,1,MPI_INT,nums.data(),1,MPI_INT,0,pgc->mpi_comm);

    size_t total=0;
    if(p->mpirank==0)
    for(int q=0; q<p->mpi_size; ++q)
    {
        counts[q] = nums[q]*NF;
        displ[q] = int(total)*NF;
        total += nums[q];
    }

    all.resize(total+1);

    MPI_Gatherv(P.data(),num*NF,MPI_DOUBLE,all.data(),counts.data(),displ.data(),MPI_DOUBLE,0,pgc->mpi_comm);

    all.resize(total);

    // model -> world coordinates
    if(p->mpirank==0 && p->cms_flag==1)
    for(auto &a : all)
    {
        double xw = p->Xout(a.x,a.y);
        double yw = p->Yout(a.x,a.y);
        double uw = a.u*cos(p->alpha_grid) - a.v*sin(p->alpha_grid);
        double vw = a.u*sin(p->alpha_grid) + a.v*cos(p->alpha_grid);
        a.x=xw; a.y=yw;
        a.u=uw; a.v=vw;
    }
}

void nhflow_particle_f::print_vtp(lexer *p, const std::vector<nhflow_particle_data> &all)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_NHFLOW_Particles/REEF3D-NHFLOW-Particles-%08i.vtp",printcount);

    std::ofstream result(name, std::ios::binary);

    const uint32_t np = uint32_t(all.size());

    // appended raw data, offsets in bytes
    const uint64_t spoint = 4 + uint64_t(np)*3*sizeof(float);
    const uint64_t sscal  = 4 + uint64_t(np)*sizeof(float);
    const uint64_t sint   = 4 + uint64_t(np)*sizeof(int32_t);

    uint64_t offset=0;
    result<<"<?xml version=\"1.0\"?>\n";
    result<<"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
    result<<"<PolyData>\n";
    result<<"<FieldData>\n<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"<<std::setprecision(12)<<p->simtime<<"</DataArray>\n</FieldData>\n";
    result<<"<Piece NumberOfPoints=\""<<np<<"\" NumberOfVerts=\""<<np<<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    result<<"<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>\n</Points>\n";
    offset+=spoint;

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=spoint;
    result<<"<DataArray type=\"Int32\" Name=\"id\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sint;
    result<<"<DataArray type=\"Int32\" Name=\"release\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sint;
    result<<"<DataArray type=\"Int32\" Name=\"state\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sint;
    result<<"<DataArray type=\"Float32\" Name=\"age\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sscal;
    result<<"</PointData>\n";

    result<<"<Verts>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sint;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sint;
    result<<"</Verts>\n";

    result<<"</Piece>\n</PolyData>\n<AppendedData encoding=\"raw\">\n_";

    auto wsize = [&](uint64_t s){ uint32_t b = uint32_t(s-4); result.write(reinterpret_cast<const char*>(&b),4); };
    auto wf = [&](double v){ float f=float(v); result.write(reinterpret_cast<const char*>(&f),sizeof(float)); };
    auto wi = [&](int32_t v){ result.write(reinterpret_cast<const char*>(&v),sizeof(int32_t)); };

    wsize(spoint);
    for(auto &a : all) { wf(a.x); wf(a.y); wf(a.z); }

    wsize(spoint);
    for(auto &a : all) { wf(a.u); wf(a.v); wf(a.w); }

    wsize(sint);
    for(auto &a : all) wi(int32_t(a.id));

    wsize(sint);
    for(auto &a : all) wi(int32_t(a.src));

    wsize(sint);
    for(auto &a : all) wi(int32_t(a.state));

    wsize(sscal);
    for(auto &a : all) wf(p->simtime-a.t0);

    wsize(sint);
    for(uint32_t q=0; q<np; ++q) wi(int32_t(q));

    wsize(sint);
    for(uint32_t q=0; q<np; ++q) wi(int32_t(q+1));

    result<<"\n</AppendedData>\n</VTKFile>\n";
    result.close();
}

void nhflow_particle_f::print_csv(lexer *p, const std::vector<nhflow_particle_data> &all)
{
    if(!csvout.is_open())
    return;

    csvout<<std::setprecision(10);
    for(auto &a : all)
    csvout<<p->simtime<<","<<(long long)(a.id)<<","<<int(a.src)<<","<<a.x<<","<<a.y<<","<<a.z<<","
          <<a.u<<","<<a.v<<","<<a.w<<","<<int(a.state)<<"\n";
    csvout.flush();
}

void nhflow_particle_f::print_log(lexer *p, ghostcell *pgc)
{
    long long cnt[4]={0,0,0,0};
    for(auto &a : P)
    ++cnt[std::min(std::max(int(a.state),0),3)];

    long long loc[5] = {cnt[0],cnt[1],cnt[2],cnt[3],numout};
    long long glob[5];
    MPI_Allreduce(loc,glob,5,MPI_LONG_LONG,MPI_SUM,pgc->mpi_comm);

    long long released=0;
    for(auto &r : R)
    released += r.released;

    if(p->mpirank==0 && logout.is_open())
    {
    logout<<std::setprecision(9)<<p->simtime<<" \t "<<released<<" \t "<<glob[0]+glob[1]+glob[2]+glob[3]<<" \t "
          <<glob[0]<<" \t "<<glob[1]<<" \t "<<glob[2]<<" \t "<<glob[3]<<" \t "<<glob[4]<<endl;
    }
}
