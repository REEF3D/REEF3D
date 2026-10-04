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

#include"fem_solid.h"
#include<fstream>
#include<ostream>
#include<iomanip>
#include<stdexcept>

void fem_solid::write_vtu(const std::string& filename) const
{
    // intact elements as hexahedra, debris particles as vertices
    std::vector<int> cells;
    for(int e=0; e<nelem(); ++e)
    if(elems[e].alive)
    cells.push_back(e);

    const int nh = (int)cells.size();
    const int no = (int)orphan.size();

    std::ofstream f(filename.c_str());
    if(!f)
    throw std::runtime_error("FEM: cannot write "+filename);

    f<<std::setprecision(7);
    f<<"<?xml version=\"1.0\"?>\n";
    f<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f<<"<UnstructuredGrid>\n";
    f<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"<<t<<"</DataArray></FieldData>\n";
    f<<"<Piece NumberOfPoints=\""<<nnode()<<"\" NumberOfCells=\""<<nh+no<<"\">\n";

    f<<"<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i=0; i<nnode(); ++i)
    f<<x[i](0)<<" "<<x[i](1)<<" "<<x[i](2)<<"\n";
    f<<"</DataArray>\n</Points>\n";

    f<<"<PointData Vectors=\"displacement\">\n";
    f<<"<DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i=0; i<nnode(); ++i)
    f<<x[i](0)-X[i](0)<<" "<<x[i](1)-X[i](1)<<" "<<x[i](2)-X[i](2)<<"\n";
    f<<"</DataArray>\n";
    f<<"<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i=0; i<nnode(); ++i)
    f<<v[i](0)<<" "<<v[i](1)<<" "<<v[i](2)<<"\n";
    f<<"</DataArray>\n";
    f<<"<DataArray type=\"Float64\" Name=\"load\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i=0; i<nnode(); ++i)
    {
        const Vec3 L = fext[i] + fcpl[i] - mfl[i]*grav;
        f<<L(0)<<" "<<L(1)<<" "<<L(2)<<"\n";
    }
    f<<"</DataArray>\n";
    f<<"<DataArray type=\"Int32\" Name=\"support\" format=\"ascii\">\n";
    for(int i=0; i<nnode(); ++i) f<<(fixed[i] ? 1 : 0)<<"\n";
    f<<"</DataArray>\n";
    f<<"</PointData>\n";

    f<<"<CellData Scalars=\"vonMises\">\n";
    f<<"<DataArray type=\"Float64\" Name=\"vonMises\" format=\"ascii\">\n";
    for(int e : cells) f<<elems[e].svm<<"\n";
    for(int q=0; q<no; ++q) f<<"0\n";
    f<<"</DataArray>\n";

    f<<"<DataArray type=\"Float64\" Name=\"damage\" format=\"ascii\">\n";
    for(int e : cells)
    {
        double d = 0.0;
        for(int g=0; g<ngp; ++g) d += gps[e*ngp+g].d;
        f<<d/double(ngp)<<"\n";
    }
    for(int q=0; q<no; ++q) f<<"1\n";
    f<<"</DataArray>\n";

    f<<"<DataArray type=\"Float64\" Name=\"plastic_strain\" format=\"ascii\">\n";
    for(int e : cells)
    {
        double p = 0.0;
        for(int g=0; g<ngp; ++g) p += gps[e*ngp+g].ep;
        f<<p/double(ngp)<<"\n";
    }
    for(int q=0; q<no; ++q) f<<"0\n";
    f<<"</DataArray>\n";

    f<<"<DataArray type=\"Float64\" Name=\"utilisation\" format=\"ascii\">\n";
    for(int e : cells) f<<elems[e].util<<"\n";
    for(int q=0; q<no; ++q) f<<"-1\n";
    f<<"</DataArray>\n";

    f<<"<DataArray type=\"Int32\" Name=\"material\" format=\"ascii\">\n";
    for(int e : cells) f<<mats[elems[e].mat].id<<"\n";
    for(int q=0; q<no; ++q) f<<"-1\n";
    f<<"</DataArray>\n";
    f<<"</CellData>\n";

    f<<"<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int e : cells)
    {
        for(int a=0; a<8; ++a) f<<elems[e].n[a]<<" ";
        f<<"\n";
    }
    for(int i : orphan) f<<i<<"\n";
    f<<"</DataArray>\n";
    f<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int o = 0;
    for(int q=0; q<nh; ++q) {o += 8; f<<o<<"\n";}
    for(int q=0; q<no; ++q) {o += 1; f<<o<<"\n";}
    f<<"</DataArray>\n";
    f<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for(int q=0; q<nh; ++q) f<<"12\n";
    for(int q=0; q<no; ++q) f<<"1\n";
    f<<"</DataArray>\n</Cells>\n";

    f<<"</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

void fem_solid::info(std::ostream& os) const
{
    os<<"FEM: "<<nelem()<<" hex8 elements ("<<(ngp==8 ? "full" : "reduced")<<" integration), "
      <<nnode()<<" nodes, lattice "<<nx<<"x"<<ny<<"x"<<nz<<" h = "<<hx<<" "<<hy<<" "<<hz
      <<", "<<faces.size()<<" surface faces, "<<nbodies<<" bodies, dt_crit "<<dtcrit<<" s"
      <<(plane_strain ? ", plane strain" : "")<<"\n";

    for(const material& mt : mats)
    {
        os<<"FEM: material "<<mt.id<<" "
          <<(mt.type==MAT_ELASTIC ? "elastic" : mt.type==MAT_J2 ? "plastic" : "concrete")
          <<" rho "<<mt.rho<<" E "<<mt.E<<" nu "<<mt.nu<<" c_p "<<mt.cp;
        if(mt.type==MAT_J2) os<<" sigma_y "<<mt.sigy<<" H "<<mt.H<<" eps_fail "<<mt.epsfail;
        if(mt.type==MAT_CONCRETE) os<<" ft "<<mt.ft<<" Gf "<<mt.Gf<<" fc "<<mt.fc<<" Gc "<<mt.Gc;
        os<<"\n";
    }

    int nfix = 0;
    for(unsigned char f : fixed) if(f) ++nfix;
    os<<"FEM: "<<nfix<<" supported nodes, mass "<<[&]{double s=0; for(double q:m) s+=q; return s;}()<<" kg\n";
}
