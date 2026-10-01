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

#include"fnpf_amr.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<mpi.h>
#include<iomanip>
#include<cstdio>

//  Output of FNPF AMR: the free surface of level 0 and of every patch as VTK rectilinear
//  grids (REEF3D_FNPF_AMR/*.vtr, indexed by a .vtm per output time, P 20 / P 30 as the FNPF
//  output), the wave gauges P 51 from the finest grid that holds them (bilinear between the
//  cell centres), and a log.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }
}

void fnpf_amr::print(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    bool doprint=false;

    if((p->count%p->P20==0 && p->P30<0.0 && p->P10>0 && p->P20>0) || (p->count==0 && p->P30<0.0))
    doprint=true;

    if((p->simtime>printtime_amr && p->P30>0.0 && p->P10>0) || (p->count==0 && p->P30>0.0))
    {
        doprint=true;
        printtime_amr += p->P30;
    }

    if(p->P51>0)
    gauges(p,c,pgc);

    if(p->mpirank==0)
    logout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<p->dt<<" \t "<<patches_total<<" \t "<<cells_total<<" \t "
          <<lap_it_last<<" \t "<<setprecision(4)<<lap_res_last<<endl;

    if(p->mpirank==0 && doprint)
    cout<<"FNPF AMR: "<<patches_total<<" patches, "<<cells_total<<" columns; time in patch stages "<<setprecision(4)<<tm[0]
        <<" s, Laplace "<<tm[1]<<" s (preconditioner "<<tm[2]<<" s, operator "<<tm[3]<<" s), mean iterations "
        <<(lap_solves>0 ? double(lap_it_total)/lap_solves : 0.0)<<endl;

    if(!doprint)
    return;

    write_vtr0(p,c);

    for(int n=0; n<(int)P.size(); ++n)
    write_vtr(p,*FP(n),n);

    int np = (int)P.size();
    vector<int> all(p->mpi_size,0);
    MPI_Allgather(&np,1,MPI_INT,&all[0],1,MPI_INT,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        char name[256];
        snprintf(name,sizeof(name),"./REEF3D_FNPF_AMR/REEF3D-FNPF-AMR-%08i.vtm",printcount_amr);
        ofstream out(name);
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n<vtkMultiBlockDataSet>\n";
        out<<"<Block index=\"0\" name=\"level 0\">\n";
        for(int r=0; r<p->mpi_size; ++r)
        out<<"<DataSet index=\""<<r<<"\" file=\"REEF3D-FNPF-AMR-L0-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<".vtr\"/>\n";
        out<<setfill(' ')<<"</Block>\n<Block index=\"1\" name=\"patches\">\n";
        int idx=0;
        for(int r=0; r<p->mpi_size; ++r)
        for(int q=0; q<all[r]; ++q)
        {
        out<<"<DataSet index=\""<<idx<<"\" file=\"REEF3D-FNPF-AMR-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<"-"<<setw(4)<<q+1<<".vtr\"/>\n";
        out<<setfill(' ');
        ++idx;
        }
        out<<"</Block>\n</vtkMultiBlockDataSet>\n</VTKFile>\n";
        out.close();
    }

    ++printcount_amr;
}

// surface elevation at the P 51 gauges from the finest grid that holds them
void fnpf_amr::gauges(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const int ng = p->P51;
    vector<double> lv(ng,-1.0), val(ng,0.0);

    auto find = [&](lexer *q, int i0, int i1, int j0, int j1, double x, double y, int &ii, int &jj)
    {
        ii=-1; jj=-1;
        for(int a=i0; a<=i1; ++a)
        if(x>=q->XN[a+marge] && x<q->XN[a+1+marge]) { ii=a; break; }
        for(int a=j0; a<=j1; ++a)
        if(y>=q->YN[a+marge] && y<q->YN[a+1+marge]) { jj=a; break; }
        return ii>=0 && jj>=0;
    };

    // bilinear between the cell centres around (x,y); (ii,jj) is the cell that holds the point
    auto ipol = [&](lexer *q, slice &f, int ii, int jj, double x, double y)
    {
        int a = (x>=q->XP[ii+marge]) ? ii : ii-1;
        int b = (y>=q->YP[jj+marge]) ? jj : jj-1;
        double wx = (x-q->XP[a+marge])/(q->XP[a+1+marge]-q->XP[a+marge]);
        double wy = (y-q->YP[b+marge])/(q->YP[b+1+marge]-q->YP[b+marge]);
        if(q->flagslice4[lij(q,a,b)]<0 || q->flagslice4[lij(q,a+1,b)]<0 || q->flagslice4[lij(q,a,b+1)]<0 || q->flagslice4[lij(q,a+1,b+1)]<0)
        return f(ii,jj);
        return (1.0-wx)*(1.0-wy)*f(a,b) + wx*(1.0-wy)*f(a+1,b) + (1.0-wx)*wy*f(a,b+1) + wx*wy*f(a+1,b+1);
    };

    for(int k=0; k<ng; ++k)
    {
        int ii,jj;
        const double x = p->P51_x[k], y = p->P51_y[k];
        if(find(p,0,NX0-1,0,NY0-1,x,y,ii,jj))
        {
            lv[k]=0.0;
            val[k]=ipol(p,c->eta,ii,jj,x,y);
        }
        for(auto q : P)
        if(q->lev>lv[k] && find(q->pp,EXT,EXT+q->nx-1,EXT,EXT+q->ny-1,x,y,ii,jj))
        {
            lv[k]=q->lev;
            val[k]=ipol(q->pp,FP(q)->c->eta,ii,jj,x,y);
        }
    }

    vector<double> lmax(ng);
    MPI_Allreduce(&lv[0],&lmax[0],ng,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    for(int k=0; k<ng; ++k)
    if(lv[k]!=lmax[k])
    val[k]=0.0;
    vector<double> vs(ng);
    MPI_Allreduce(&val[0],&vs[0],ng,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        if(!gaugeout.is_open())
        {
            gaugeout.open("./REEF3D_FNPF_AMR/REEF3D_FNPF_AMR_gauges.dat");
            gaugeout<<"# simtime";
            for(int k=0; k<ng; ++k)
            gaugeout<<" \t eta("<<p->P51_x[k]<<","<<p->P51_y[k]<<")";
            gaugeout<<" \t levels"<<endl;
        }
        gaugeout<<setprecision(10)<<p->simtime;
        for(int k=0; k<ng; ++k)
        gaugeout<<" \t "<<setprecision(10)<<vs[k];
        for(int k=0; k<ng; ++k)
        gaugeout<<" \t "<<int(lmax[k]);
        gaugeout<<endl;
    }
}

void fnpf_amr::write_vtr0(lexer *p, fdm_fnpf *c)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_FNPF_AMR/REEF3D-FNPF-AMR-L0-%08i-%04i.vtr",printcount_amr,p->mpirank+1);

    const int nx=p->knox, ny=p->knoy, m=marge;
    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<RectilinearGrid WholeExtent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">0</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=0; jj<ny; ++jj)
        {
            for(int ii=0; ii<nx; ++ii)
            out<<setprecision(17)<<f(ii,jj)+shift<<" ";
            out<<"\n";
        }
        out<<"</DataArray>\n";
    };
    field("eta",c->eta,0.0);
    field("elevation",c->eta,p->wd);
    field("Fifsf",c->Fifsf,0.0);
    field("Fz",c->Fz,0.0);
    field("bed",c->bed,0.0);

    out<<"</CellData>\n<Coordinates>\n";
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=0; ii<=nx; ++ii) out<<setprecision(12)<<p->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=0; jj<=ny; ++jj) out<<setprecision(12)<<p->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}

void fnpf_amr::write_vtr(lexer *p, fnpf_amr_patch &c, int id)
{
    lexer *pp = c.pp;
    fdm_fnpf *cc = c.c;
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_FNPF_AMR/REEF3D-FNPF-AMR-%08i-%04i-%04i.vtr",printcount_amr,p->mpirank+1,id+1);

    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<RectilinearGrid WholeExtent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<p->simtime<<"</DataArray>";
    out<<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\">"<<c.lev<<"</DataArray></FieldData>\n";
    out<<"<Piece Extent=\"0 "<<c.nx<<" 0 "<<c.ny<<" 0 0\">\n";
    out<<"<CellData Scalars=\"eta\">\n";

    auto field = [&](const char *nm, slice &f, double shift)
    {
        out<<"<DataArray type=\"Float64\" Name=\""<<nm<<"\" format=\"ascii\">\n";
        for(int jj=EXT; jj<EXT+c.ny; ++jj)
        {
            for(int ii=EXT; ii<EXT+c.nx; ++ii)
            out<<setprecision(17)<<f(ii,jj)+shift<<" ";
            out<<"\n";
        }
        out<<"</DataArray>\n";
    };
    field("eta",cc->eta,0.0);
    field("elevation",cc->eta,p->wd);
    field("Fifsf",cc->Fifsf,0.0);
    field("Fz",cc->Fz,0.0);
    field("bed",cc->bed,0.0);

    out<<"</CellData>\n<Coordinates>\n";
    const int m = marge;
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=EXT; ii<=EXT+c.nx; ++ii) out<<setprecision(12)<<pp->XN[ii+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=EXT; jj<=EXT+c.ny; ++jj) out<<setprecision(12)<<pp->YN[jj+m]<<" ";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">0</DataArray>\n";
    out<<"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
    out.close();
}
