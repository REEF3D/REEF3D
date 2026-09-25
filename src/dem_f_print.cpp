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

void dem_f::print(lexer *p, ghostcell *pgc)
{
    print_log(p,pgc);

    if(p->E15<=0.0 || p->simtime<printtime-1.0e-12)
    return;

    if(p->mpirank==0)
    {
        print_vtp(p);
        print_state(p);
    }

    ++printcount;
    printtime += p->E15;
}

void dem_f::print_log(lexer *p, ghostcell *pgc)
{
    if(p->mpirank!=0 || p->count%std::max(1,p->P12)!=0)
    return;

    int nact=0;
    for(auto &B : core.bodies)
    nact += B.active ? 1 : 0;

    std::streamsize prec = cout.precision();
    cout<<"DEM: particles "<<nact<<" substeps "<<nsub<<" contacts "<<core.ncontacts<<" iterations "<<maxiter_used
        <<" residual "<<setprecision(3)<<core.residual<<" max pen "<<core.maxpen<<" KE "<<core.kinetic_energy()
        <<" time "<<steptime<<setprecision(prec)<<endl;
}

void dem_f::print_vtp(lexer *p)
{
    char name[200];
    snprintf(name,sizeof(name),"./REEF3D_DEM_VTP/REEF3D-DEM-%06i.vtp",printcount);

    ofstream out(name);
    if(!out.is_open())
    return;

    int np=0, nt=0;
    for(auto &B : core.bodies)
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
    for(auto &B : core.bodies)
    if(B.active)
    for(auto &v : core.shapes[B.shape].vert)
    {
        dem_vec x = B.x + B.R*v;
        out<<x(0)<<" "<<x(1)<<" "<<x(2)<<"\n";
    }
    out<<"</DataArray>\n</Points>\n";

    out<<"<PointData>\n";
    out<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(auto &B : core.bodies)
    if(B.active)
    for(auto &v : core.shapes[B.shape].vert)
    {
        dem_vec u = B.v + B.w.cross(B.R*v);
        out<<u(0)<<" "<<u(1)<<" "<<u(2)<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"id\" format=\"ascii\">\n";
    for(auto &B : core.bodies)
    if(B.active)
    for(size_t q=0; q<core.shapes[B.shape].vert.size(); ++q)
    out<<B.id<<"\n";
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"resolved\" format=\"ascii\">\n";
    for(auto &B : core.bodies)
    if(B.active)
    for(size_t q=0; q<core.shapes[B.shape].vert.size(); ++q)
    out<<B.mode<<"\n";
    out<<"</DataArray>\n";
    out<<"</PointData>\n";

    out<<"<Polys>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int offset=0;
    for(auto &B : core.bodies)
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
}

void dem_f::print_state(lexer *p)
{
    ofstream out("./REEF3D_DEM/REEF3D-DEM-state.dat", printcount==0 ? ios::out : ios::app);
    if(!out.is_open())
    return;

    if(printcount==0)
    out<<"# time id x y z qw qx qy qz u v w wx wy wz active resolved Fx Fy Fz ufx ufy ufz"<<endl;

    out<<setprecision(10);
    for(auto &B : core.bodies)
    out<<p->simtime<<" "<<B.id<<" "<<B.x(0)<<" "<<B.x(1)<<" "<<B.x(2)<<" "
       <<B.q.w()<<" "<<B.q.x()<<" "<<B.q.y()<<" "<<B.q.z()<<" "
       <<B.v(0)<<" "<<B.v(1)<<" "<<B.v(2)<<" "<<B.w(0)<<" "<<B.w(1)<<" "<<B.w(2)<<" "
       <<(B.active?1:0)<<" "<<B.mode<<" "
       <<Fhyd[B.id](0)<<" "<<Fhyd[B.id](1)<<" "<<Fhyd[B.id](2)<<" "
       <<ufl[B.id](0)<<" "<<ufl[B.id](1)<<" "<<ufl[B.id](2)<<"\n";
}
