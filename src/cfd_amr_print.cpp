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

#include"cfd_amr.h"
#include"reef3d_mpmd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"lagoon_output.h"
#include<cstdio>
#include<vector>
#include<sys/stat.h>
#include<sys/types.h>

//  Output of the patches (P 10 > 0): at the output times of printer_CFD (P 20, P 30 and the time
//  windows P 35), one VTK rectilinear grid
//  (.vtr) per patch in REEF3D_CFD_AMR/, and a multiblock file (.vtm) with all patches of the output
//  (open it next to the level-0 output).  Cell data: level set, pressure, velocity (cell centre), the
//  level of the patch and a flag for the cells under a finer patch.
//  P 18: the patches go into the LAGOON store too, as the AMR set cfd_amr of 3D blocks (rank 0
//  gathers them, lagoon_amr_output); P 18 1: instead of the .vtr and .vtm files.

namespace
{
void write_floats(FILE *fp, const std::vector<float> &v)
{
    int n=0;
    for(float x : v)
    {
        fprintf(fp,"%.7g",x);
        fputc((++n%8==0) ? '\n' : ' ',fp);
    }
    fputc('\n',fp);
}
}

void cfd_amr::print(lexer *p, ghostcell *pgc, bool after)
{
    if(p->P10==0)
    return;

    const bool by_count = (p->P20>0 && p->P30<0.0 && p->P34<0.0 && p->count%p->P20==0);
    const bool by_time = (p->P30>0.0 && p->P34<0.0 && (p->simtime>p->printtime || p->count==0));
    if(after && p->count!=0)
    return;

    // P 35: the time windows ts te dt, the same test as printer_CFD (called before it in the step,
    // so the output number is that of its level-0 output).  The next output time of a window starts
    // at ts, after a restart (I 45) at the first window time not before the start.
    bool by_window = false;
    if(p->P35>0 && !after)
    {
        if((int)p35_next.size()!=p->P35)
        {
            p35_next.assign(p->P35,0.0);
            for(int qn=0; qn<p->P35; ++qn)
            {
                p35_next[qn] = p->P35_ts[qn];
                if(p->P35_dt[qn]>0.0)
                while(p35_next[qn]<p->simtime)
                p35_next[qn] += p->P35_dt[qn];
            }
        }

        for(int qn=0; qn<p->P35; ++qn)
        if(p->simtime>p35_next[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
        {
            by_window = true;
            p35_next[qn] += p->P35_dt[qn];
        }
    }

    if(!by_count && !by_time && !by_window)
    return;

    // after: called after printer_CFD at the start (its output counter has moved on)
    const int num = (p->P15==2) ? p->count : p->printcount - (after ? 1 : 0);

    // P 18: the patches in the LAGOON store; P 18 1: instead of the .vtr and .vtm files
    bool stored = false;
    if(p->P18>0)
    stored = print_lagoon(p,pgc,num);
    if(!lagoon_amr_output::files_needed(p,stored))
    return;

    if(myrank==0)
    mkdir("./REEF3D_CFD_AMR",0777);
    MPI_Barrier(reef3d_comm());

    char name[256];

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        fdm *a = c->a;
        const int nx=pp->knox, ny=pp->knoy, nz=pp->knoz;
        const size_t nc = (size_t)nx*ny*nz;

        snprintf(name,sizeof(name),"./REEF3D_CFD_AMR/REEF3D-CFD-AMR-%08i-%05i.vtr",num,c->gid);
        FILE *fp = fopen(name,"w");
        if(fp==nullptr)
        continue;

        fprintf(fp,"<?xml version=\"1.0\"?>\n<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fp,"<RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n<Piece Extent=\"0 %d 0 %d 0 %d\">\n",nx,ny,nz,nx,ny,nz);
        fprintf(fp,"<FieldData>\n<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\"> %.10g </DataArray>\n</FieldData>\n",p->simtime);
        fprintf(fp,"<CellData Scalars=\"phi\" Vectors=\"velocity\">\n");

        std::vector<float> v(nc), vel(3*nc);
        auto put = [&](const char *nm, auto fn)
        {
            size_t m=0;
            for(int k=0; k<nz; ++k)
            for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
            v[m++] = (float)fn(i,j,k);
            fprintf(fp,"<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n",nm);
            write_floats(fp,v);
            fprintf(fp,"</DataArray>\n");
        };

        put("phi",[&](int i, int j, int k) { return a->phi(i,j,k); });
        put("pressure",[&](int i, int j, int k) { return a->press(i,j,k); });
        put("level",[&](int i, int j, int k) { return double(c->lev); });
        put("covered",[&](int i, int j, int k) { const int I[3] = {i+c->lo[0],j+c->lo[1],k+c->lo[2]}; return covered(c->lev,I) ? 1.0 : 0.0; });
        if(scal)
        {
            field &S = sfield(id);
            put(scal==1 ? "T" : "C",[&](int i, int j, int k) { return S(i,j,k); });
        }

        size_t m=0;
        for(int k=0; k<nz; ++k)
        for(int j=0; j<ny; ++j)
        for(int i=0; i<nx; ++i)
        {
            vel[m++] = (float)(0.5*(a->u(i,j,k)+a->u(i-1,j,k)));
            vel[m++] = (float)(0.5*(a->v(i,j,k)+a->v(i,j-1,k)));
            vel[m++] = (float)(0.5*(a->w(i,j,k)+a->w(i,j,k-1)));
        }
        fprintf(fp,"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        write_floats(fp,vel);
        fprintf(fp,"</DataArray>\n</CellData>\n<Coordinates>\n");

        const double *xn[3] = {pp->XN,pp->YN,pp->ZN};
        const int n3[3] = {nx,ny,nz};
        for(int d=0; d<3; ++d)
        {
            std::vector<float> x(n3[d]+1);
            for(int q=0; q<=n3[d]; ++q)
            x[q] = (float)xn[d][q+marge];
            fprintf(fp,"<DataArray type=\"Float32\" format=\"ascii\">\n");
            write_floats(fp,x);
            fprintf(fp,"</DataArray>\n");
        }
        fprintf(fp,"</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
        fclose(fp);
    }

    if(myrank==0)
    {
        snprintf(name,sizeof(name),"./REEF3D_CFD_AMR/REEF3D-CFD-AMR-%08i.vtm",num);
        FILE *fp = fopen(name,"w");
        if(fp!=nullptr)
        {
            fprintf(fp,"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n<vtkMultiBlockDataSet>\n");
            for(size_t g=0; g<GP.size(); ++g)
            fprintf(fp,"<DataSet index=\"%d\" file=\"REEF3D-CFD-AMR-%08i-%05i.vtr\"/>\n",(int)g,num,(int)g);
            fprintf(fp,"</vtkMultiBlockDataSet>\n</VTKFile>\n");
            fclose(fp);
        }
    }
}

// The patches of this rank as 3D blocks for the LAGOON store (cell values as the .vtr files have
// them: x fastest, then y, then z; the velocity at the cell centres as u, v, w), gathered on rank 0
// and written to the AMR set cfd_amr of ./REEF3D_CFD.lagoon. True on every rank when the output is
// in the store.
bool cfd_amr::print_lagoon(lexer *p, ghostcell *pgc, int num)
{
    static lagoon_amr_output *writer = nullptr;

    std::vector<lagoon_amr::field> fields = {{"phi",false}, {"pressure",false}, {"level",true}, {"covered",true},
                                             {"u",false}, {"v",false}, {"w",false}};
    if(scal)
    fields.push_back({scal==1 ? "T" : "C",false});

    std::vector<lagoon_amr::grid> grids;
    for(int id=0; id<(int)P.size(); ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;
        fdm *a = c->a;
        lagoon_amr::grid g;
        g.level = c->lev;
        g.nx = pp->knox;
        g.ny = pp->knoy;
        g.nz = pp->knoz;
        const double *xn[3] = {pp->XN,pp->YN,pp->ZN};
        std::vector<double> *axes[3] = {&g.x,&g.y,&g.z};
        const int n3[3] = {g.nx,g.ny,g.nz};
        for(int d=0; d<3; ++d)
        for(int q=0; q<=n3[d]; ++q)
        axes[d]->push_back(xn[d][q+marge]);

        const size_t nc = (size_t)g.nx*g.ny*g.nz;
        g.values.reserve(fields.size()*nc);
        auto put = [&](auto fn)
        {
            for(int k=0; k<g.nz; ++k)
            for(int j=0; j<g.ny; ++j)
            for(int i=0; i<g.nx; ++i)
            g.values.push_back(fn(i,j,k));
        };
        put([&](int i, int j, int k) { return a->phi(i,j,k); });
        put([&](int i, int j, int k) { return a->press(i,j,k); });
        put([&](int i, int j, int k) { return double(c->lev); });
        put([&](int i, int j, int k) { const int I[3] = {i+c->lo[0],j+c->lo[1],k+c->lo[2]}; return covered(c->lev,I) ? 1.0 : 0.0; });
        put([&](int i, int j, int k) { return 0.5*(a->u(i,j,k)+a->u(i-1,j,k)); });
        put([&](int i, int j, int k) { return 0.5*(a->v(i,j,k)+a->v(i,j-1,k)); });
        put([&](int i, int j, int k) { return 0.5*(a->w(i,j,k)+a->w(i,j,k-1)); });
        if(scal)
        {
            field &S = sfield(id);
            put([&](int i, int j, int k) { return S(i,j,k); });
        }
        grids.push_back(std::move(g));
    }

    if(writer==nullptr)
    writer = new lagoon_amr_output("CFD",fields);
    return writer->write(p,pgc,grids,num);
}
