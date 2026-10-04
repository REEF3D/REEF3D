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

#include"dem_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include<fstream>
#include<iomanip>
#include<cstdio>
#include<mpi.h>
#include"runlog.h"
#include"lagoon_store.h"

namespace
{
// the run's LAGOON store and run record (P 18)
std::string lagoon_solver(lexer *p)
{
    return p->A10==6 ? "CFD" : p->A10==3 ? "FNPF" : p->A10==2 ? "SFLOW" : "NHFLOW";
}

std::string lagoon_run(lexer *p)
{
    return p->plog ? "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}" : "";
}
}

void dem_f::print(lexer *p, ghostcell *pgc)
{
    print_log(p,pgc);

    if(p->E15<=0.0 || p->simtime<printtime-1.0e-12)
    return;

    // collective: owned particles are gathered on rank 0
    vector<dem_body> all;
    gather_output(pgc,all);

    if(p->mpirank==0)
    {
        // P 18: the elements in the LAGOON store as well (dem_dem: points, velocity,
        // id, resolved and the triangles, as in the VTP file); P 18 2: instead of it
        bool stored = false;
        if(p->P18>0)
        {
            static lagoon_particles *writer = nullptr;
            if(writer==nullptr)
            {
                const std::string solver = lagoon_solver(p);
                writer = new lagoon_particles("./REEF3D_" + solver + ".lagoon", solver, "dem_dem", "dem",
                                              "REEF3D_DEM_VTP", lagoon_run(p),
                                              {{"velocity", 3, false}, {"id", 1, true}, {"resolved", 1, true}},
                                              1, "polys");
            }
            std::vector<float> xyz, velocity;
            std::vector<int32_t> id, resolved, connectivity, offsets;
            int offset = 0;
            for(auto &B : all)
            if(B.active)
            {
                const dem_shape &S = core.shapes[B.shape];
                for(auto &v : S.vert)
                {
                    const dem_vec x = B.x + B.R*v;
                    const dem_vec u = B.v + B.w.cross(B.R*v);
                    for(int c=0; c<3; ++c)
                    {
                        xyz.push_back(float(x(c)));
                        velocity.push_back(float(u(c)));
                    }
                    id.push_back(int32_t(B.id));
                    resolved.push_back(int32_t(B.mode));
                }
                for(size_t t=0; t<S.tri.size(); t+=3)
                {
                    for(int k=0; k<3; ++k)
                    connectivity.push_back(int32_t(S.tri[t+k] + offset));
                    offsets.push_back(int32_t(connectivity.size()));
                }
                offset += int(S.vert.size());
            }
            stored = writer->output(p->simtime, printcount, id.size(), xyz.data(),
                                    {velocity.data(), id.data(), resolved.data()}, connectivity, offsets);
        }
        if(!(stored && p->P18==2))
        print_vtp(p,all);
        print_state(p,all);
    }

    ++printcount;
    printtime += p->E15;
}

void dem_f::print_log(lexer *p, ghostcell *pgc)
{
    if(p->count%std::max(1,p->P12)!=0)
    return;

    // collective statistics
    double sum[3] = {0.0,0.0,0.0};     // active particles, contacts, kinetic energy
    double mx[4] = {double(maxiter_used),core.residual,core.maxpen,steptime};

    for(auto &B : core.bodies)
    if(!B.ghost && B.active && (B.tier==0 || p->mpirank==0))
    {
        sum[0] += 1.0;
        if(!B.fixed)
        {
            dem_vec wb = B.R.transpose()*B.w;
            sum[2] += 0.5*B.m*B.v.squaredNorm() + 0.5*wb.dot(B.Ib.cwiseProduct(wb));
        }
    }
    sum[1] = core.ncontacts;

    MPI_Allreduce(MPI_IN_PLACE,sum,3,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,mx,4,MPI_DOUBLE,MPI_MAX,pgc->mpi_comm);

    if(p->mpirank!=0)
    return;

    std::streamsize prec = cout.precision();
    cout<<"DEM: particles "<<int(sum[0])<<" substeps "<<nsub<<" contacts "<<int(sum[1])<<" iterations "<<int(mx[0])
        <<" residual "<<setprecision(3)<<mx[1]<<" max pen "<<mx[2]<<" KE "<<sum[2]
        <<" time "<<mx[3]<<setprecision(prec)<<endl;
}

void dem_f::print_vtp(lexer *p, const vector<dem_body> &bodies)
{
    char name[200];
    snprintf(name,sizeof(name),"./REEF3D_DEM_VTP/REEF3D-DEM-%06i.vtp",printcount);

    ofstream out(name);
    if(!out.is_open())
    return;

    int np=0, nt=0;
    for(auto &B : bodies)
    if(B.active)
    {
        np += core.shapes[B.shape].vert.size();
        nt += core.shapes[B.shape].tri.size()/3;
    }

    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";
    out<<"<FieldData>\n<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<setprecision(12)<<p->simtime<<"</DataArray>\n</FieldData>\n";
    out<<"<Piece NumberOfPoints=\""<<np<<"\" NumberOfPolys=\""<<nt<<"\">\n";

    out<<setprecision(8);
    out<<"<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(auto &B : bodies)
    if(B.active)
    for(auto &v : core.shapes[B.shape].vert)
    {
        dem_vec x = B.x + B.R*v;
        out<<x(0)<<" "<<x(1)<<" "<<x(2)<<"\n";
    }
    out<<"</DataArray>\n</Points>\n";

    out<<"<PointData>\n";
    out<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(auto &B : bodies)
    if(B.active)
    for(auto &v : core.shapes[B.shape].vert)
    {
        dem_vec u = B.v + B.w.cross(B.R*v);
        out<<u(0)<<" "<<u(1)<<" "<<u(2)<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"id\" format=\"ascii\">\n";
    for(auto &B : bodies)
    if(B.active)
    for(size_t q=0; q<core.shapes[B.shape].vert.size(); ++q)
    out<<B.id<<"\n";
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"resolved\" format=\"ascii\">\n";
    for(auto &B : bodies)
    if(B.active)
    for(size_t q=0; q<core.shapes[B.shape].vert.size(); ++q)
    out<<B.mode<<"\n";
    out<<"</DataArray>\n";
    out<<"</PointData>\n";

    out<<"<Polys>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int offset=0;
    for(auto &B : bodies)
    if(B.active)
    {
        const dem_shape &S = core.shapes[B.shape];
        for(size_t t=0; t<S.tri.size(); t+=3)
        out<<S.tri[t]+offset<<" "<<S.tri[t+1]+offset<<" "<<S.tri[t+2]+offset<<"\n";
        offset += S.vert.size();
    }
    out<<"</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(int t=1; t<=nt; ++t)
    out<<3*t<<"\n";
    out<<"</DataArray>\n</Polys>\n";

    out<<"</Piece>\n</PolyData>\n</VTKFile>\n";
    out.close();

    if(p->plog)
    p->plog->written(p,printcount,"dem","dem",name,0);
}

void dem_f::print_state(lexer *p, const vector<dem_body> &bodies)
{
    ofstream out("./REEF3D_DEM/REEF3D-DEM-state.dat", printcount==0 ? ios::out : ios::app);
    if(!out.is_open())
    return;

    if(printcount==0)
    out<<"# time id x y z qw qx qy qz u v w wx wy wz active resolved Fx Fy Fz ufx ufy ufz"<<endl;

    out<<setprecision(10);
    for(auto &B : bodies)
    out<<p->simtime<<" "<<B.id<<" "<<B.x(0)<<" "<<B.x(1)<<" "<<B.x(2)<<" "
       <<B.q.w()<<" "<<B.q.x()<<" "<<B.q.y()<<" "<<B.q.z()<<" "
       <<B.v(0)<<" "<<B.v(1)<<" "<<B.v(2)<<" "<<B.w(0)<<" "<<B.w(1)<<" "<<B.w(2)<<" "
       <<(B.active?1:0)<<" "<<B.mode<<" "
       <<B.cpl.Fhyd(0)<<" "<<B.cpl.Fhyd(1)<<" "<<B.cpl.Fhyd(2)<<" "
       <<B.cpl.ufl(0)<<" "<<B.cpl.ufl(1)<<" "<<B.cpl.ufl(2)<<"\n";
}
