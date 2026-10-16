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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cstdio>
#include<vector>
#include<sys/stat.h>
#include<sys/types.h>

//  Output of the patches (P 10 > 0): at the output times of printer_CFD, one VTK rectilinear grid
//  (.vtr) per patch in REEF3D_CFD_AMR/, and a multiblock file (.vtm) with all patches of the output
//  (open it next to the level-0 output).  Cell data: level set, pressure, velocity (cell centre), the
//  level of the patch and a flag for the cells under a finer patch.

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

void cfd_amr::print(lexer *p, bool after)
{
    if(p->P10==0)
    return;

    const bool by_count = (p->P20>0 && p->P30<0.0 && p->P34<0.0 && p->count%p->P20==0);
    const bool by_time = (p->P30>0.0 && p->P34<0.0 && (p->simtime>p->printtime || p->count==0));
    if(after && p->count!=0)
    return;
    if(!by_count && !by_time)
    return;

    // after: called after printer_CFD at the start (its output counter has moved on)
    const int num = (p->P15==2) ? p->count : p->printcount - (after ? 1 : 0);

    if(myrank==0)
    mkdir("./REEF3D_CFD_AMR",0777);
    MPI_Barrier(MPI_COMM_WORLD);

    char name[256];

    for(int id=0; id<(int)P.size(); ++id)
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
