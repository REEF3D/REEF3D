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

#include"fnpf_vtp_fsf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<sstream>
#include<vector>
#include"lagoon_output.h"

namespace
{
lagoon_surface *lagoon_free_surface = nullptr;  // P 18
}

fnpf_vtp_fsf::fnpf_vtp_fsf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->I40==0)
        p->printtime=0.0;

    printcount=0;

    // Create Folder
    if(p->mpirank==0)
        mkdir("./REEF3D_FNPF_VTP_FSF",0777);

    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    // 2D
    if(p->j_dir==0)
    {
        gcval_eta = 155;
        gcval_fifsf = 160;
    }
}

void fnpf_vtp_fsf::start(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{
    print2D(p,c,pgc);
}

void fnpf_vtp_fsf::print2D(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    pgc->gcsl_start4(p,c->test2D,1);

    SLICELOOP4
    {
        if(c->breaking(i,j)>=1)
        c->breaking_print(i,j)=double(c->breaking(i,j));
        
        else if(c->breaking(i,j)==0)
        c->breaking_print(i,j)=0.0;
    }

    pgc->gcsl_start4(p,c->breaking_print,50);

    int num=0;
    if(p->P15==1)
    num = printcount;
    
    else if(p->P15==2)
    num = p->count;

    if(p->mpirank==0 && lagoon_output::vtu_files(p))  // P 18 1: the store instead
    pvtp(p,num);

    // offsets
    n=0;
    offset[n]=0;
    ++n;

    // Points
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;

    // velocity
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;

    // Fifsf
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // elevation
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // depth
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // breaking
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // coastline
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // test
    if(p->P23==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
        ++n;
    }

    // Hs
    if(p->P110==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
        ++n;
    }

    // Cells
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    // the piece in memory (lagoon_surface::piece writes the file where it is wanted); every array
    // made in a buffer and written at once: the same bytes as value by value
    snprintf(name,sizeof(name),"./REEF3D_FNPF_VTP_FSF/REEF3D-FNPF-FSF-%08i-%06i.vtp",num,p->mpirank+1);
    std::ostringstream result(std::ios::out | std::ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result, offset, n);

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"Fifsf\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"eta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"depth\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"coastline\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P23==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    if(p->P110==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"Hs\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    result<<"</PointData>\n";

    vtp3D::polys(result, offset, n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    std::vector<float> fb;
    fb.reserve(size_t(p->pointnum2D)*3);
    auto put_size = [&](int bytes){ result.write((char*)&bytes, sizeof(int)); };
    auto put_floats = [&](){ result.write((char*)fb.data(), std::streamsize(sizeof(float)*fb.size())); fb.clear(); };
    // a point array of the slice f, interpolated to the nodes
    auto put_slice = [&](slice &f)
    {
        put_size(sizeof(float)*p->pointnum2D);
        TPSLICELOOP
        fb.push_back(float(p->sl_ipol4(f)));
        put_floats();
    };

    //  XYZ
    put_size(sizeof(float)*p->pointnum2D*3);
    TPSLICELOOP
    {
        fb.push_back(float(p->Xout(p->XN[IP1],p->YN[JP1])));
        fb.push_back(float(p->Yout(p->XN[IP1],p->YN[JP1])));
        fb.push_back(float(p->nhf_ipol4eta(p->wet,c->eta,c->bed)+p->wd));
    }
    put_floats();

    //  Velocities
    put_size(sizeof(float)*p->pointnum2D*3);
    TPSLICELOOP
    {
        k = p->knoz-1;
        
        // without the 3D velocities: the surface ones (fnpf_memory.h), the same node
        if(c->U==nullptr)
        {
        fb.push_back(float(c->Us(i,j)));
        fb.push_back(float(c->Vs(i,j)));
        fb.push_back(float(c->Ws(i,j)));
        continue;
        }

        float ffn=float(c->U[FIJKp1]);
        if(k==-1 && j==-1)
            ffn=float(c->U[FIJp1Kp1]);
        fb.push_back(ffn);

        ffn=float(c->V[FIJKp1]);
        if(k==-1 && j==-1)
            ffn=float(c->V[FIJp1Kp1]);
        fb.push_back(ffn);

        ffn=float(c->W[FIJKp1]);
        if(k==-1 && j==-1)
            ffn=float(c->W[FIJp1Kp1]);
        fb.push_back(ffn);
    }
    put_floats();

    put_slice(c->Fifsf);
    put_slice(c->eta);
    put_slice(c->depth);
    put_slice(c->breaking_print);
    put_slice(c->coastline);

    if(p->P23==1)
    put_slice(c->test2D);

    if(p->P110==1)
    put_slice(c->Hs);

    //  Connectivity
    std::vector<int> ib;
    ib.reserve(size_t(p->polygon_sum)*3);
    put_size(sizeof(int)*p->polygon_sum*3);
    SLICEBASELOOP
    {
        // Triangle 1
        ib.push_back(int(c->nodeval2D(i-1,j-1))-1);
        ib.push_back(int(c->nodeval2D(i,j-1))-1);
        ib.push_back(int(c->nodeval2D(i,j))-1);

        // Triangle 2
        ib.push_back(int(c->nodeval2D(i-1,j-1))-1);
        ib.push_back(int(c->nodeval2D(i,j))-1);
        ib.push_back(int(c->nodeval2D(i-1,j))-1);
    }
    result.write((char*)ib.data(), std::streamsize(sizeof(int)*ib.size()));

    //  Offset of Connectivity
    ib.clear();
    put_size(sizeof(int)*p->polygon_sum);
    for(n=0;n<p->polygon_sum;++n)
    ib.push_back((n+1)*3);
    result.write((char*)ib.data(), std::streamsize(sizeof(int)*ib.size()));

    vtp3D::footer(result);

    // the file (P 18 0, 2) and the piece in the LAGOON store (P 18 1: instead of the file, P 18 2: as well)
    std::string piece = std::move(result).str();
    lagoon_surface::piece(p,pgc,lagoon_free_surface,"FNPF","free_surface","REEF3D_FNPF_VTP_FSF",name,num,piece);

    ++printcount;
}
