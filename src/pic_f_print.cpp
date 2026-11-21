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

#include"pic_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"lagoon_store.h"
#include"runlog.h"
#include<mpi.h>
#include<cstdio>
#include<cstdint>
#include<iomanip>
#include<iostream>

// Output of the particles, all gathered on rank 0:
// P 18 1: the LAGOON store ./REEF3D_CFD.lagoon, particle set cfd_pic (the VTP files only if
//         the store cannot be written), P 18 2: the store and the VTP files, P 18 0: the VTP
//         files ./REEF3D_CFD_PIC/REEF3D-CFD-PIC-<output>.vtp
// Point fields: velocity, id (int32), rank (the MPI rank that has the particle).
// When: every F 87 seconds; F 87 <= 0: with the CFD volume output (P 30 seconds or every
// P 20 iterations). REEF3D_CFD_PIC/REEF3D-CFD-PIC-Log.dat: the number of particles, those
// that left the domain, those sent between the ranks, those added and removed in the water
// (pic_f::reseed) and the time of the free surface from the particles [s], at every output.

double pic_f::print_interval(lexer *p) const
{
    if(p->F87>0.0)
    return p->F87;

    if(p->P30>0.0)
    return p->P30;

    return -1.0;
}

bool pic_f::print_due(lexer *p) const
{
    const double eps = 1.0e-10;
    const double dtp = print_interval(p);

    if(dtp>0.0)
    return p->simtime >= printtime - eps;

    if(p->P20>0)
    return p->count%p->P20==0;

    return false;
}

void pic_f::print(lexer *p, ghostcell *pgc)
{
    const double eps = 1.0e-10;

    std::vector<double> all;
    gather(p,pgc,all);

    if(p->mpirank==0)
    {
        bool stored = false;
        if(p->P18>0)
        stored = print_lagoon(p,all);

        if(!(stored && p->P18==1))
        print_vtp(p,all);

        const long long n = (long long)(all.size()/NG);

        logout<<setprecision(9)<<p->simtime<<" \t "<<n<<" \t "<<lost_total<<" \t "<<sent_total
              <<" \t "<<added_total<<" \t "<<removed_total<<" \t "<<phitime<<endl;

        cout<<"PIC output "<<printcount<<": "<<n<<" particles, "<<lost_total<<" left the domain"
            <<(stored ? " (LAGOON store)" : "")<<endl;
    }

    ++printcount;

    const double dtp = print_interval(p);
    if(dtp>0.0)
    while(printtime <= p->simtime + eps)
    printtime += dtp;
}

// rank 0: x y z u v w id rank of every particle
void pic_f::gather(lexer *p, ghostcell *pgc, std::vector<double> &all)
{
    std::vector<double> mine;
    mine.reserve(P.size()*NG);

    for(const auto &q : P)
    {
        const double row[NG] = {q.x,q.y,q.z, q.u,q.v,q.w, q.id, double(p->mpirank)};
        mine.insert(mine.end(),row,row+NG);
    }

    int num = int(mine.size());
    std::vector<int> nums(p->mpi_size,0), displ(p->mpi_size,0);

    MPI_Gather(&num,1,MPI_INT,nums.data(),1,MPI_INT,0,pgc->mpi_comm);

    size_t total=0;
    if(p->mpirank==0)
    for(int r=0; r<p->mpi_size; ++r)
    {
        displ[r] = int(total);
        total += size_t(nums[r]);
    }

    all.resize(total+1);
    MPI_Gatherv(mine.data(),num,MPI_DOUBLE,all.data(),nums.data(),displ.data(),MPI_DOUBLE,0,pgc->mpi_comm);
    all.resize(total);
}

// rank 0: one output of the particle set cfd_pic; false if it could not be written
bool pic_f::print_lagoon(lexer *p, const std::vector<double> &all)
{
    static lagoon_particles *writer = nullptr;

    if(writer==nullptr)
    {
        std::string run;
        if(p->plog)
        run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";

        const std::vector<lagoon_particles::field> fields = {{"velocity",3,false}, {"id",1,true}, {"rank",1,true}};
        writer = new lagoon_particles("./REEF3D_CFD.lagoon", "CFD", "cfd_pic", "particles",
                                      "REEF3D_CFD_PIC", run, fields);
    }

    if(!writer->usable())
    return false;

    const size_t n = all.size()/NG;
    std::vector<float> xyz(3*n), vel(3*n);
    std::vector<int32_t> id(n), rank(n);

    for(size_t q=0; q<n; ++q)
    {
        const double *r = &all[q*NG];
        for(int c=0; c<3; ++c)
        {
            xyz[3*q+c] = float(r[c]);
            vel[3*q+c] = float(r[3+c]);
        }
        id[q] = int32_t(int64_t(r[6]));
        rank[q] = int32_t(r[7]);
    }

    return writer->output(p->simtime, printcount, n, xyz.data(), {vel.data(), id.data(), rank.data()});
}

// rank 0: one VTP file of all particles
void pic_f::print_vtp(lexer *p, const std::vector<double> &all)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_CFD_PIC/REEF3D-CFD-PIC-%08i.vtp",printcount);

    std::ofstream result(name, std::ios::binary);

    const uint32_t np = uint32_t(all.size()/NG);

    // appended raw data, offsets in bytes
    const uint64_t spoint = 4 + uint64_t(np)*3*sizeof(float);
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
    result<<"<DataArray type=\"Int32\" Name=\"rank\" format=\"appended\" offset=\""<<offset<<"\"/>\n";
    offset+=sint;
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
    for(uint32_t q=0; q<np; ++q) { wf(all[q*NG]); wf(all[q*NG+1]); wf(all[q*NG+2]); }

    wsize(spoint);
    for(uint32_t q=0; q<np; ++q) { wf(all[q*NG+3]); wf(all[q*NG+4]); wf(all[q*NG+5]); }

    wsize(sint);
    for(uint32_t q=0; q<np; ++q) wi(int32_t(int64_t(all[q*NG+6])));

    wsize(sint);
    for(uint32_t q=0; q<np; ++q) wi(int32_t(all[q*NG+7]));

    wsize(sint);
    for(uint32_t q=0; q<np; ++q) wi(int32_t(q));

    wsize(sint);
    for(uint32_t q=0; q<np; ++q) wi(int32_t(q+1));

    result<<"\n</AppendedData>\n</VTKFile>\n";
    result.close();
}
