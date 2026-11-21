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

#include"soil_stress.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<cstdio>
#include<cstdint>
#include<iostream>
#include<iomanip>

using namespace std;

// cell values of the soil: displacement (3), seepage velocity (3, Darcy), excess pore pressure,
// effective stresses sxx syy szz txy tyz txz, liquefaction potential
// LP = (p - p_b) / ((1+2 K0)/3 gamma' d), d the depth below the bed, p_b the pore pressure at the bed
// of the column (Jeng 1997, Zen and Yamazaki 1990): LP >= 1 momentarily liquefied
void soil_stress::postprocess(lexer *p)
{
    const vget g{X.data(),-1,0};
    std::fill(post.begin(),post.end(),0.0);

    // bed of the columns: the highest face soil-fluid, its pore pressure (over the ranks of the column)
    const int ncol = knox*knoy;
    std::vector<double> zbed(ncol,-1.0e30), pbed(ncol,0.0);

    for(int i=0;i<knox;++i)
    for(int j=0;j<knoy;++j)
    for(int k=0;k<knoz;++k)
    if(type(i,j,k)==SOIL && type(i,j,k+1)==FREE)
    {
        zbed[i*knoy+j] = p->ZN[k+1+p->marge];
        pbed[i*knoy+j] = pb[id(i,j,k+1)];
    }

    {
        std::vector<double> zmax(ncol);
        MPI_Allreduce(zbed.data(),zmax.data(),ncol,MPI_DOUBLE,MPI_MAX,colcomm);
        for(int n=0;n<ncol;++n)
        if(zbed[n]!=zmax[n])
        pbed[n] = 0.0;
        std::vector<double> psum(ncol);
        MPI_Allreduce(pbed.data(),psum.data(),ncol,MPI_DOUBLE,MPI_SUM,colcomm);
        zbed = zmax;
        pbed = psum;
    }

    auto flux = [&](int i, int j, int k, int dir) -> double
    {
        // Darcy flux through the face dir (1 x+, 3 y+, 5 z+; 0, 2, 4 the minus faces) of the soil cell, out of it
        int d = dir;
        const double T = flow_trans(i,j,k,d);
        if(T==0.0)
        return 0.0;
        int ii=i, jj=j, kk=k;
        double area;
        if(dir<2) {ii += dir==0?-1:1; area = DYN(j)*DZN(k);}
        else if(dir<4) {jj += dir==2?-1:1; area = DXN(i)*DZN(k);}
        else {kk += dir==4?-1:1; area = DXN(i)*DYN(j);}
        const double pn = d==SOIL ? X[3*N+id(ii,jj,kk)] : pb[id(ii,jj,kk)];
        return T*(X[3*N+id(i,j,k)] - pn)/area;
    };

    const double lpfac = (1.0+2.0*K0)/3.0*gsub;

    for(int i=0;i<knox;++i)
    for(int j=0;j<knoy;++j)
    for(int k=0;k<knoz;++k)
    if(type(i,j,k)==SOIL)
    {
        const int q = id(i,j,k);
        auto out = [&](int m, double val) {post[size_t(m)*N+q] = val;};

        out(0,0.5*(U(g,0,i,j,k)+U(g,0,i-1,j,k)));
        out(1,dim3 ? 0.5*(U(g,1,i,j,k)+U(g,1,i,j-1,k)) : 0.0);
        out(2,0.5*(U(g,2,i,j,k)+U(g,2,i,j,k-1)));

        out(3,0.5*(flux(i,j,k,1) - flux(i,j,k,0)));
        out(4,dim3 ? 0.5*(flux(i,j,k,3) - flux(i,j,k,2)) : 0.0);
        out(5,0.5*(flux(i,j,k,5) - flux(i,j,k,4)));

        const double pc = X[3*N+q];
        out(6,pc);

        double sxx,syy,szz,d1,d2;
        normal_stress(i,j,k,g,sxx,syy,szz);
        out(7,sxx);
        out(8,syy);
        out(9,szz);

        if(dim3)
        {
        out(10,0.25*(edge_xy(i-1,j-1,k,g,d1,d2)+edge_xy(i,j-1,k,g,d1,d2)+edge_xy(i-1,j,k,g,d1,d2)+edge_xy(i,j,k,g,d1,d2)));
        out(11,0.25*(edge_yz(i,j-1,k-1,g,d1,d2)+edge_yz(i,j,k-1,g,d1,d2)+edge_yz(i,j-1,k,g,d1,d2)+edge_yz(i,j,k,g,d1,d2)));
        }
        out(12,0.25*(edge_xz(i-1,j,k-1,g,d1,d2)+edge_xz(i,j,k-1,g,d1,d2)+edge_xz(i-1,j,k,g,d1,d2)+edge_xz(i,j,k,g,d1,d2)));

        const int n = i*knoy+j;
        const double depth = zbed[n] - p->ZP[k+p->marge];
        if(zbed[n]>-1.0e29 && depth>0.0 && lpfac>0.0)
        out(13,(pc - pbed[n])/(lpfac*depth));
    }
}

void soil_stress::print_probes(lexer *p, ghostcell *pgc)
{
    const int np = p->S216;
    std::vector<double> val(size_t(5)*np,0.0), cnt(np,0.0);

    for(int n=0;n<np;++n)
    {
        int ci=-1, cj=0, ck=-1;
        for(int i=0;i<knox;++i)
        if(p->S216_x[n]>=p->XN[i+p->marge] && p->S216_x[n]<p->XN[i+1+p->marge])
        ci = i;
        if(dim3)
        {
        cj = -1;
        for(int j=0;j<knoy;++j)
        if(p->S216_y[n]>=p->YN[j+p->marge] && p->S216_y[n]<p->YN[j+1+p->marge])
        cj = j;
        }
        for(int k=0;k<knoz;++k)
        if(p->S216_z[n]>=p->ZN[k+p->marge] && p->S216_z[n]<p->ZN[k+1+p->marge])
        ck = k;

        if(ci<0 || cj<0 || ck<0)
        continue;

        const int q = id(ci,cj,ck);
        cnt[n] = 1.0;
        if(type(ci,cj,ck)!=SOIL)
        continue;
        val[5*n+0] = post[size_t(6)*N+q];
        val[5*n+1] = post[size_t(0)*N+q];
        val[5*n+2] = post[size_t(1)*N+q];
        val[5*n+3] = post[size_t(2)*N+q];
        val[5*n+4] = post[size_t(13)*N+q];
    }

    std::vector<double> vsum(val.size()), csum(np);
    MPI_Allreduce(val.data(),vsum.data(),int(val.size()),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(cnt.data(),csum.data(),np,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    if(p->mpirank==0 && probefile.is_open())
    {
        probefile<<setprecision(9)<<p->simtime;
        for(int n=0;n<np;++n)
        for(int m=0;m<5;++m)
        probefile<<"\t"<<setprecision(9)<<(csum[n]>0.0 ? vsum[5*n+m]/csum[n] : 0.0);
        probefile<<endl;
    }
}

// the soil cells as hexahedra, cell data; REEF3D_CFD_Soil/REEF3D-CFD-Soil-<n>.pvtu with a piece per rank
void soil_stress::print_vtu(lexer *p)
{
    std::vector<int> cells;
    for(int i=0;i<knox;++i)
    for(int j=0;j<knoy;++j)
    for(int k=0;k<knoz;++k)
    if(type(i,j,k)==SOIL)
    cells.push_back(id(i,j,k));

    const int nc = int(cells.size());
    const int npts = 8*nc;

    static const char *names[] = {"displacement","seepage_velocity","pore_pressure","sigma_xx","sigma_yy","sigma_zz",
                                  "tau_xy","tau_yz","tau_xz","liquefaction_potential"};
    static const int ncomp[] = {3,3,1,1,1,1,1,1,1,1};
    static const int first[] = {0,3,6,7,8,9,10,11,12,13};
    const int narr = 10;

    char name[300];
    snprintf(name,sizeof(name),"./REEF3D_CFD_Soil/REEF3D-CFD-Soil-%08i-%06i.vtu",printcount,p->mpirank+1);
    ofstream result(name,ios::binary);

    // offsets of the appended data: points, connectivity, offsets, types, cell arrays
    std::vector<uint64_t> off;
    uint64_t o = 0;
    auto add = [&](uint64_t bytes) {off.push_back(o); o += sizeof(uint32_t) + bytes;};
    add(uint64_t(npts)*3*sizeof(float));
    add(uint64_t(npts)*sizeof(int32_t));
    add(uint64_t(nc)*sizeof(int32_t));
    add(uint64_t(nc)*sizeof(uint8_t));
    for(int m=0;m<narr;++m)
    add(uint64_t(nc)*ncomp[m]*sizeof(float));

    result<<"<?xml version=\"1.0\"?>"<<endl;
    result<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">"<<endl;
    result<<"<UnstructuredGrid>"<<endl;
    result<<"<FieldData>"<<endl<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<setprecision(12)<<p->simtime<<endl<<"</DataArray>"<<endl<<"</FieldData>"<<endl;
    result<<"<Piece NumberOfPoints=\""<<npts<<"\" NumberOfCells=\""<<nc<<"\">"<<endl;
    result<<"<Points>"<<endl<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<off[0]<<"\"/>"<<endl<<"</Points>"<<endl;
    result<<"<Cells>"<<endl;
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<off[1]<<"\"/>"<<endl;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<off[2]<<"\"/>"<<endl;
    result<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<<off[3]<<"\"/>"<<endl;
    result<<"</Cells>"<<endl;
    result<<"<CellData>"<<endl;
    for(int m=0;m<narr;++m)
    result<<"<DataArray type=\"Float32\" Name=\""<<names[m]<<"\" NumberOfComponents=\""<<ncomp[m]<<"\" format=\"appended\" offset=\""<<off[4+m]<<"\"/>"<<endl;
    result<<"</CellData>"<<endl;
    result<<"</Piece>"<<endl<<"</UnstructuredGrid>"<<endl;
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

    auto block = [&](uint64_t bytes) {uint32_t h = uint32_t(bytes); result.write((char*)&h,sizeof(h));};
    float ff;
    int32_t ii;

    block(uint64_t(npts)*3*sizeof(float));
    for(int q : cells)
    {
        const int i = q/(NY*NZ) - 1;
        const int j = (q/NZ)%NY - 1;
        const int k = q%NZ - 1;
        const double xs[2] = {p->XN[i+p->marge],p->XN[i+1+p->marge]};
        const double ys[2] = {p->YN[j+p->marge],p->YN[j+1+p->marge]};
        const double zs[2] = {p->ZN[k+p->marge],p->ZN[k+1+p->marge]};
        static const int cx[8] = {0,1,1,0,0,1,1,0};
        static const int cy[8] = {0,0,1,1,0,0,1,1};
        static const int cz[8] = {0,0,0,0,1,1,1,1};
        for(int c=0;c<8;++c)
        {
            ff = float(xs[cx[c]]); result.write((char*)&ff,sizeof(float));
            ff = float(ys[cy[c]]); result.write((char*)&ff,sizeof(float));
            ff = float(zs[cz[c]]); result.write((char*)&ff,sizeof(float));
        }
    }

    block(uint64_t(npts)*sizeof(int32_t));
    for(int n=0;n<npts;++n)
    {
        ii = n;
        result.write((char*)&ii,sizeof(int32_t));
    }

    block(uint64_t(nc)*sizeof(int32_t));
    for(int n=0;n<nc;++n)
    {
        ii = 8*(n+1);
        result.write((char*)&ii,sizeof(int32_t));
    }

    block(uint64_t(nc)*sizeof(uint8_t));
    for(int n=0;n<nc;++n)
    {
        const uint8_t t = 12;
        result.write((char*)&t,sizeof(uint8_t));
    }

    for(int m=0;m<narr;++m)
    {
        block(uint64_t(nc)*ncomp[m]*sizeof(float));
        for(int q : cells)
        for(int c=0;c<ncomp[m];++c)
        {
            ff = float(post[size_t(first[m]+c)*N+q]);
            result.write((char*)&ff,sizeof(float));
        }
    }
    result<<endl<<"</AppendedData>"<<endl<<"</VTKFile>"<<endl;
    result.close();

    // the parallel file
    if(p->mpirank==0)
    {
        snprintf(name,sizeof(name),"./REEF3D_CFD_Soil/REEF3D-CFD-Soil-%08i.pvtu",printcount);
        ofstream pv(name);
        pv<<"<?xml version=\"1.0\"?>"<<endl;
        pv<<"<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">"<<endl;
        pv<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;
        pv<<"<PPoints>"<<endl<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl<<"</PPoints>"<<endl;
        pv<<"<PCellData>"<<endl;
        for(int m=0;m<narr;++m)
        pv<<"<PDataArray type=\"Float32\" Name=\""<<names[m]<<"\" NumberOfComponents=\""<<ncomp[m]<<"\"/>"<<endl;
        pv<<"</PCellData>"<<endl;
        for(int r=0;r<p->mpi_size;++r)
        {
            snprintf(name,sizeof(name),"REEF3D-CFD-Soil-%08i-%06i.vtu",printcount,r+1);
            pv<<"<Piece Source=\""<<name<<"\"/>"<<endl;
        }
        pv<<"</PUnstructuredGrid>"<<endl<<"</VTKFile>"<<endl;
    }

    ++printcount;
}
