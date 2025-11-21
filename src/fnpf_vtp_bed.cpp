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

#include"fnpf_vtp_bed.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_vtp_bed::fnpf_vtp_bed(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->I40==0)
        p->printtime=0.0;

    printcount=0;

    // Create Folder
    if(p->mpirank==0)
        mkdir("./REEF3D_FNPF_VTP_BED",0777);
}

void fnpf_vtp_bed::start(lexer *p, fdm_fnpf *c, ghostcell* pgc, ioflow *pflow)
{
    print2D(p,c,pgc);
}

void fnpf_vtp_bed::print2D(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{
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

    // elevation
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // depth
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // Cells
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    // Open File
    ofstream result;
    sprintf(name,"./REEF3D_FNPF_VTP_BED/REEF3D-FNPF-BED-%08i-%06i.vtp",num,p->mpirank+1);
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"depth\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";

    vtp3D::polys(result,offset,n);

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

        ffn=float(p->sl_ipol4(c->bed));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Elevation
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(c->bed));
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

    //  Connectivity
    iin=sizeof(int)*(p->polygon_sum)*3;
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
    iin=sizeof(int)*(p->polygon_sum);
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
