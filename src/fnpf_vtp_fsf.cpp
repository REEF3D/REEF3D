/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

    if(p->mpirank==0)
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

    // Open File
    ofstream result;
    sprintf(name,"./REEF3D_FNPF_VTP_FSF/REEF3D-FNPF-FSF-%08i-%06i.vtp",num,p->mpirank+1);
    result.open(name, ios::binary);

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

    float ffn;
    int iin;
    //  XYZ
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->XN[IP1]);
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->YN[JP1]);
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->nhf_ipol4eta(p->wet,c->eta,c->bed)+p->wd);
        result.write((char*)&ffn, sizeof(float));
    }

    //  Velocities
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        k = p->knoz-1;

        ffn=float(c->U[FIJKp1]);
        if(k==-1 && j==-1)
            ffn=float(c->U[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof(float));

        ffn=float(c->V[FIJKp1]);
        if(k==-1 && j==-1)
            ffn=float(c->V[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof(float));

        ffn=float(c->W[FIJKp1]);
        if(k==-1 && j==-1)
            ffn=float(c->W[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof(float));
    }

    //  Fifsf
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(c->Fifsf));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Eta
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(c->eta));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Depth
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(c->depth));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Breaking
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(c->breaking_print));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Coastline
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(c->coastline));
        result.write((char*)&ffn, sizeof(float));
    }

    //  test
    if(p->P23==1)
    {
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(c->test2D));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  Hs
    if(p->P110==1)
    {
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(c->Hs));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  Connectivity
    iin=sizeof(int)*p->polygon_sum*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
        // Triangle 1
        iin=int(c->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(c->nodeval2D(i,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(c->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        // Triangle 2
        iin=int(c->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(c->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(c->nodeval2D(i-1,j))-1;
        result.write((char*)&iin, sizeof(int));
    }

    //  Offset of Connectivity
    iin=sizeof(int)*p->polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<p->polygon_sum;++n)
    {
        iin=(n+1)*3;
        result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);

    result.close();

    ++printcount;
}
