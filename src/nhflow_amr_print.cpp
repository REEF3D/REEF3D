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

#include"nhflow_amr.h"
#include"nhflow_amr_fill.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"lagoon_output.h"
#include"vtr3D.h"
#include<mpi.h>
#include<iomanip>
#include<cstdio>

//  Output of NHFLOW AMR: the free surface of level 0 and of every patch as VTK rectilinear grids
//  (REEF3D_NHFLOW_AMR/*.vtr, indexed by a .vtm per output time, P 20 / P 30 as the NHFLOW
//  output), the wave gauges P 51 from the finest grid that holds them (bilinear between the cell
//  centres), and a log with the water volume of the leaf cells. P 18: the grids also in the
//  LAGOON store (P 18 1: instead of the .vtr and .vtm files).

using namespace nhflow_amr_detail;

// water volume of the leaf cells (level-0 and patch cells not covered by a finer patch)
double nhflow_amr::mass(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double v=0.0;

    for(int ii=0; ii<NX0; ++ii)
    for(int jj=0; jj<NY0; ++jj)
    {
        if(p->flagslice4[lij(p,ii,jj)]<0)
        continue;
        if(maxlev>=1 && covered(1,2*(ii+O0i),2*(jj+O0j)))
        continue;
        v += d->WL(ii,jj)*p->DXN[ii+marge]*p->DYN[jj+marge];
    }

    for(auto q : P)
    {
        lexer *pp = q->pp;
        fdm_nhf *dd = NP(q)->d;
        for(int ii=EXT; ii<EXT+q->nx; ++ii)
        for(int jj=EXT; jj<EXT+q->ny; ++jj)
        {
            if(pp->flagslice4[lij(pp,ii,jj)]<0)
            continue;
            const int I = ii-EXT+q->I0, J = jj-EXT+q->J0;
            if(q->lev<maxlev && covered(q->lev+1,2*I,2*J))
            continue;
            v += dd->WL(ii,jj)*pp->DXN[ii+marge]*pp->DYN[jj+marge];
        }
    }

    return pgc->globalsum(v);
}

void nhflow_amr::print(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    bool doprint=false;

    if((p->count%p->P20==0 && p->P30<0.0 && p->P10>0 && p->P20>0) || (p->count==0 && p->P30<0.0))
    doprint=true;

    if((p->simtime>printtime_amr && p->P30>0.0) || (p->count==0 && p->P30>0.0))
    {
        doprint=true;
        printtime_amr += p->P30;
    }

    if(p->P51>0)
    gauges(p,d,pgc);

    const double m = mass(p,d,pgc);

    if(p->mpirank==0)
    logout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<p->dt<<" \t "<<patches_total<<" \t "<<cells_total<<" \t "
          <<pr_it_last<<" \t "<<setprecision(4)<<pr_res_last<<" \t "<<setprecision(15)<<m<<" \t "<<setprecision(6)<<(m0>0.0 ? (m-m0)/m0 : 0.0)<<" \t "<<layout_id<<endl;

    if(p->mpirank==0 && (doprint || p->count%500==0))
    cout<<"NHFLOW AMR: "<<patches_total<<" patches, "<<cells_total<<" columns; time fill "<<setprecision(4)<<tm[0]
        <<" s, patch stages "<<tm[1]<<" s, flux exchange "<<tm[2]<<" s, restriction "<<tm[3]<<" s, pressure "<<tm[4]
        <<" s (preconditioner "<<tm[5]<<" s, operator "<<tm[6]<<" s), mean iterations "
        <<(pr_solves>0 ? double(pr_it_total)/pr_solves : 0.0);
    if(p->mpirank==0 && (doprint || p->count%500==0) && sub==1)
    cout<<"; G 7 1: level solves "<<sub_lv_n<<" (mean iterations "<<setprecision(3)<<(sub_lv_n>0 ? double(sub_lv_it)/sub_lv_n : 0.0)
        <<"), synchronisation "<<sub_sy_n<<" (mean iterations "<<(sub_sy_n>0 ? double(sub_sy_it)/sub_sy_n : 0.0)<<", "<<setprecision(4)<<tsync<<" s)";
    if(p->mpirank==0 && (doprint || p->count%500==0) && cdiff)
    cout<<", diffusion "<<setprecision(4)<<tdiff<<" s, mean iterations "<<(df_solves>0 ? double(df_it_total)/df_solves : 0.0);
    if(p->mpirank==0 && (doprint || p->count%500==0) && regrid_int>0)
    cout<<", regrid "<<setprecision(4)<<tm[7]<<" s, layouts "<<layout_id<<", regrids skipped "<<regrids_skipped;
    // load per rank (several ranks): the cells computed on the rank (level 0 including its covered
    // columns, the patches without their EXT cells) and the time in the patch stages, max over mean
    // of the ranks - what balancing the patches over the ranks could gain at most
    if((doprint || p->count%500==0) && p->mpi_size>1)
    {
        double cl = double(p0->knox)*p0->knoy*p0->knoz;
        for(auto q : P)
        cl += double(q->nx)*q->ny*q->pp->knoz;

        double v[2] = {cl, tm[1]};
        double vmax[2] = {cl, tm[1]};
        pgc->globalmax(vmax,2);
        const double cmean = pgc->globalsum(v[0])/p->mpi_size;
        const double tmean = pgc->globalsum(v[1])/p->mpi_size;

        if(p->mpirank==0)
        cout<<"; load max/mean over "<<p->mpi_size<<" ranks: cells "<<setprecision(3)<<(cmean>0.0 ? vmax[0]/cmean : 1.0)
            <<", patch stages "<<(tmean>0.0 ? vmax[1]/tmean : 1.0);
    }

    if(p->mpirank==0 && (doprint || p->count%500==0))
    cout<<"; water volume change "<<setprecision(3)<<(m0>0.0 ? (m-m0)/m0 : 0.0)<<endl;

    if(!doprint)
    return;

    // P 18: the grids in the LAGOON store; P 18 1: instead of the .vtr and .vtm files
    bool stored = false;
    if(p->P18>0)
    stored = print_lagoon(p,d,pgc);

    if(!lagoon_amr_output::files_needed(p,stored))
    {
        ++printcount_amr;
        return;
    }

    write_vtr0(p,d);

    for(int n=0; n<(int)P.size(); ++n)
    write_vtr(p,*NP(n),n);

    // multiblock index: one block per level, level 0 of every rank + the patches of each level
    int np = (int)P.size();
    vector<int> all(p->mpi_size,0);
    MPI_Allgather(&np,1,MPI_INT,&all[0],1,MPI_INT,MPI_COMM_WORLD);

    vector<int> mylev(np), off(p->mpi_size,0);
    for(int n=0; n<np; ++n)
    mylev[n] = P[n]->lev;
    for(int r=1; r<p->mpi_size; ++r)
    off[r] = off[r-1]+all[r-1];
    vector<int> plev(off[p->mpi_size-1]+all[p->mpi_size-1]+1);
    MPI_Gatherv(np>0?&mylev[0]:NULL,np,MPI_INT,&plev[0],&all[0],&off[0],MPI_INT,0,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        int ltop=0;
        for(int r=0; r<p->mpi_size; ++r)
        for(int q=0; q<all[r]; ++q)
        ltop = max(ltop,plev[off[r]+q]);

        char name[256];
        snprintf(name,sizeof(name),"./REEF3D_NHFLOW_AMR/REEF3D-NHFLOW-AMR-%08i.vtm",printcount_amr);
        ofstream out(name);
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n<vtkMultiBlockDataSet>\n";
        out<<"<Block index=\"0\" name=\"level 0\">\n";
        for(int r=0; r<p->mpi_size; ++r)
        out<<"<DataSet index=\""<<r<<"\" file=\"REEF3D-NHFLOW-AMR-L0-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<".vtr\"/>\n";
        out<<setfill(' ')<<"</Block>\n";
        for(int l=1; l<=ltop; ++l)
        {
            out<<"<Block index=\""<<l<<"\" name=\"level "<<l<<"\">\n";
            int idx=0;
            for(int r=0; r<p->mpi_size; ++r)
            for(int q=0; q<all[r]; ++q)
            if(plev[off[r]+q]==l)
            {
                out<<"<DataSet index=\""<<idx<<"\" file=\"REEF3D-NHFLOW-AMR-"<<setw(8)<<setfill('0')<<printcount_amr<<"-"<<setw(4)<<r+1<<"-"<<setw(4)<<q+1<<".vtr\"/>\n";
                out<<setfill(' ');
                ++idx;
            }
            out<<"</Block>\n";
        }
        out<<"</vtkMultiBlockDataSet>\n</VTKFile>\n";
        out.close();
    }

    ++printcount_amr;
}

// surface elevation at the P 51 gauges from the finest grid that holds them
void nhflow_amr::gauges(lexer *p, fdm_nhf *d, ghostcell *pgc)
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
            val[k]=ipol(p,d->eta,ii,jj,x,y);
        }
        for(auto q : P)
        if(q->lev>lv[k] && find(q->pp,EXT,EXT+q->nx-1,EXT,EXT+q->ny-1,x,y,ii,jj))
        {
            lv[k]=q->lev;
            val[k]=ipol(q->pp,NP(q)->d->eta,ii,jj,x,y);
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
            gaugeout.open("./REEF3D_NHFLOW_AMR/REEF3D_NHFLOW_AMR_gauges.dat");
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

namespace
{
// the cell fields of a grid: eta, elevation, bed, and the velocities of the top layer;
// field(name, value(ii,jj)) for each, in the order of the .vtr files (and the LAGOON store)
template<class FIELD>
void amr_fields(lexer *q, fdm_nhf *d, bool wet, FIELD field)
{
    const int kt = q->knoz-1;
    auto c3 = [&](int ii, int jj) { return (ii-q->imin)*q->jmax*q->kmax + (jj-q->jmin)*q->kmax + kt-q->kmin; };

    field("eta",[&](int ii, int jj) { return d->eta(ii,jj); });
    field("elevation",[&](int ii, int jj) { return d->eta(ii,jj)+q->wd; });
    field("bed",[&](int ii, int jj) { return d->bed(ii,jj); });
    field("u_top",[&](int ii, int jj) { return d->U[c3(ii,jj)]; });
    field("v_top",[&](int ii, int jj) { return d->V[c3(ii,jj)]; });
    // G 30: the wet flag and the water depth
    if(wet)
    {
    field("wet",[&](int ii, int jj) { return double(q->wet[(ii-q->imin)*q->jmax + (jj-q->jmin)]); });
    field("WL",[&](int ii, int jj) { return d->WL(ii,jj); });
    }
}
}

void nhflow_amr::write_vtr0(lexer *p, fdm_nhf *d)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_NHFLOW_AMR/REEF3D-NHFLOW-AMR-L0-%08i-%04i.vtr",printcount_amr,p->mpirank+1);

    write_vtr_grid(p,p,d,name,0,p->knox,0,p->knoy,0);
}

void nhflow_amr::write_vtr(lexer *p, nhflow_amr_patch &c, int id)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_NHFLOW_AMR/REEF3D-NHFLOW-AMR-%08i-%04i-%04i.vtr",printcount_amr,p->mpirank+1,id+1);

    write_vtr_grid(p,c.pp,c.d,name,EXT,c.nx,EXT,c.ny,c.lev);
}

// cells [i0,i0+nx) x [j0,j0+ny) of a grid (lexer q, fields f), binary appended Float32 via vtr3D:
// the fields of amr_fields
void nhflow_amr::write_vtr_grid(lexer *p, lexer *q, fdm_nhf *f, const char *name, int i0, int nx, int j0, int ny, int lev)
{
    vtr3D vtr;
    const int m = marge;
    const int ncell = nx*ny;
    const double z0 = 0.0;

    vector<const char*> fname;
    vector<float> values;
    amr_fields(q,f,shore,[&](const char *nm, auto g)
    {
        fname.push_back(nm);
        for(int jj=j0; jj<j0+ny; ++jj)
        for(int ii=i0; ii<i0+nx; ++ii)
        values.push_back(float(g(ii,jj)));
    });
    const int nf = fname.size();

    // fields, then the coordinates
    vector<int> offset(nf+4);
    int n=0;
    offset[n]=0;
    ++n;
    for(int r=0; r<nf; ++r)
    {
        offset[n]=offset[n-1]+sizeof(float)*ncell+sizeof(int);
        ++n;
    }
    vtr.offset(offset.data(),n,nx+1,ny+1,1);

    stringstream result;
    const int ext[6] = {0,nx,0,ny,0,0};
    vtr.beginning(result,ext,p->simtime,lev);
    n=0;
    result<<"<CellData Scalars=\"eta\">\n";
    for(int r=0; r<nf; ++r)
    {
        result<<"<DataArray type=\"Float32\" Name=\""<<fname[r]<<"\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    result<<"</CellData>\n";
    vtr.ending(result,offset.data(),n);

    size_t pos = result.str().length();
    const size_t total = pos + offset[n] + 27;
    vector<char> buffer(total);
    memcpy(&buffer[0],result.str().data(),pos);

    const int iin = sizeof(float)*ncell;
    for(int r=0; r<nf; ++r)
    {
        memcpy(&buffer[pos],&iin,sizeof(int));
        pos+=sizeof(int);
        memcpy(&buffer[pos],&values[size_t(r)*ncell],sizeof(float)*ncell);
        pos+=sizeof(float)*ncell;
    }

    vtr.structureWrite(buffer,pos,&q->XN[i0+m],nx+1,&q->YN[j0+m],ny+1,&z0,1);

    FILE *file = fopen(name,"wb");
    if(file)
    {
        fwrite(buffer.data(),buffer.size(),1,file);
        fclose(file);
    }
}

// P 18: this rank's grids as its .vtr files have them (amr_fields), gathered into the LAGOON
// store; true when the output is there
bool nhflow_amr::print_lagoon(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    static lagoon_amr_output *writer = nullptr;
    const int m = marge;
    vector<lagoon_amr::grid> grids;
    vector<lagoon_amr::field> fields;
    auto take = [&](lexer *q, fdm_nhf *dd, int i0, int j0, int nx, int ny, int level)
    {
        lagoon_amr::grid g;
        g.level = level;
        g.nx = nx;
        g.ny = ny;
        for(int ii=i0; ii<=i0+nx; ++ii) g.x.push_back(q->XN[ii+m]);
        for(int jj=j0; jj<=j0+ny; ++jj) g.y.push_back(q->YN[jj+m]);
        fields.clear();
        amr_fields(q,dd,shore,[&](const char *nm, auto f)
        {
            fields.push_back({nm,false});
            for(int jj=j0; jj<j0+ny; ++jj)
            for(int ii=i0; ii<i0+nx; ++ii)
            g.values.push_back(f(ii,jj));
        });
        grids.push_back(std::move(g));
    };
    take(p,d,0,0,p->knox,p->knoy,0);
    for(int n=0; n<(int)P.size(); ++n)
    take(NP(n)->pp,NP(n)->d,EXT,EXT,NP(n)->nx,NP(n)->ny,NP(n)->lev);

    if(writer==nullptr)
    writer = new lagoon_amr_output("NHFLOW",fields);
    return writer->write(p,pgc,grids,printcount_amr);
}
