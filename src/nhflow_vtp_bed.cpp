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

#include"nhflow_vtp_bed.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_vtp_bed::nhflow_vtp_bed(lexer *p)
{
    if(p->I40==0)
        p->printtime=0.0;

    printcount=0;

    // Create Folder
    if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_VTP_BED",0777);
}

void nhflow_vtp_bed::start(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    print2D(p,d,pgc,psed);
}

void nhflow_vtp_bed::print2D(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    int num=0;
    if(p->P15==1)
        num = printcount;
    else if(p->P15==2)
        num = p->count;

    if(p->mpirank==0)
        pvtp(p,psed,num);

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
    // sediment bedlaod
    if(p->P76==1)
        psed->offset_ParaView_2D_bedload(p,offset,n);
    // sediment parameters 1
    if(p->P77==1)
        psed->offset_ParaView_2D_parameter1(p,offset,n);
    // sediment parameters 2
    if(p->P78==1)
        psed->offset_ParaView_2D_parameter2(p,offset,n);
    // bed shear stress
    if(p->P79>=1)
        psed->offset_ParaView_2D_bedshear(p,offset,n);
    // test
    if(p->P23==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
        ++n;
    }
    // Cells
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    sprintf(name,"./REEF3D_NHFLOW_VTP_BED/REEF3D-NHFLOW-BED-%08i-%06i.vtp",num,p->mpirank+1);

    // Open File
    ofstream result;
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"depth\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P76==1)
        psed->name_ParaView_bedload(p,result,offset,n);
    if(p->P77==1)
        psed->name_ParaView_parameter1(p,result,offset,n);
    if(p->P78==1)
        psed->name_ParaView_parameter2(p,result,offset,n);
    if(p->P79>=1)
        psed->name_ParaView_bedshear(p,result,offset,n);
    if(p->P23==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    result<<"</PointData>\n";

    vtp3D::polys(result,offset,n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    //  XYZ
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->XN[IP1]);
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->YN[JP1]);
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->sl_ipol4(d->bed));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Elevation
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(d->bed));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Depth
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(d->depth));
        result.write((char*)&ffn, sizeof(float));
    }

    //  sediment bedload
    if(p->P76==1)
        psed->print_2D_bedload(p,pgc,result);

    //  sediment parameter 1
    if(p->P77==1)
        psed->print_2D_parameter1(p,pgc,result);

    //  sediment parameter 2
    if(p->P78==1)
        psed->print_2D_parameter2(p,pgc,result);

    //  bed shear stress
    if(p->P79>=1)
        psed->print_2D_bedshear(p,pgc,result);

    //  Test
    if(p->P23==1)
    {
        pgc->gcsl_start4(p,d->test2D,1);

        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(d->test2D));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  kin and eps
    //pnhfturb->print_2D(p,d,pgc,result,1);

    //  Connectivity
    iin=sizeof(int)*p->polygon_sum*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
        // Triangle 1
        iin=int(d->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));


        // Triangle 2
        iin=int(d->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i-1,j))-1;
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
