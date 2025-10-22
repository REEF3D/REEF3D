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

#include"sflow_vtp_bed.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_vtp_bed::sflow_vtp_bed(lexer *p)
{
    if(p->I40==0)
        printbedtime=0.0;

    printbedcount=0;

    // Create Folder
    if(p->mpirank==0)
        mkdir("./REEF3D_SFLOW_VTP_BED",0777);
}

void sflow_vtp_bed::start(lexer *p, fdm2D* b, ghostcell* pgc, sediment *psed)
{
    // Print out based on iteration
    if((p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)  || (p->count==0 &&  p->P30<0.0))
    {
        print2D(p,b,pgc,psed);
    }

    // Print out based on time
    if((p->simtime>printbedtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
    {
        print2D(p,b,pgc,psed);

        printbedtime+=p->P30;
    }
}

void sflow_vtp_bed::print2D(lexer *p, fdm2D* b, ghostcell* pgc, sediment *psed)
{
    int num=0;
    if(p->P15==1)
        num = printbedcount;
    if(p->P15==2)
        num = p->count;

    if(p->mpirank==0)
        pvtp(p,psed,num);

    pgc->gcsl_start4(p,b->depth,50);
    pgc->gcsl_start4(p,b->bed,50);
    pgc->gcsl_start4(p,b->test,50);

    // bednode upate
    TPSLICELOOP
    {
        b->bednode(i,j) = 0.25*(b->bed(i,j) + b->bed(i+1,j) + b->bed(i,j+1) + b->bed(i,j));

        if(p->flagslice4[Im1Jm1]<0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]>0)
            b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j) + b->bed(i,j+1) + b->bed(i,j));

        if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]<0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]>0)
            b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i,j+1) + b->bed(i,j));

        if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]<0 && p->flagslice4[IJ]>0)
            b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i+1,j) + b->bed(i,j));

        if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]<0)
            b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i+1,j) + b->bed(i,j+1));
    }

    i=-1;
    j=-1;
    b->bednode(i,j) = b->bed(i+1,j+1);

    i=p->knox-1;
    j=-1;
    b->bednode(i,j) = b->bed(i,j+1);

    i=p->knox-1;
    j=p->knoy-1;
    b->bednode(i,j) = b->bed(i,j);

    i=-1;
    j=p->knoy-1;
    b->bednode(i,j) = b->bed(i+1,j);

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

    // pressure
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    // elevation
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
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
        ++n;
    }

    // Cells
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    // Open File
    sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-%08i-%06i.vtp",num,p->mpirank+1);
    ofstream result;
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P76==1)
        psed->name_ParaView_bedload(p,pgc,result,offset,n);
    if(p->P77==1)
        psed->name_ParaView_parameter1(p,pgc,result,offset,n);
    if(p->P78==1)
        psed->name_ParaView_parameter2(p,pgc,result,offset,n);
    if(p->P79>=1)
        psed->name_ParaView_bedshear(p,pgc,result,offset,n);
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
        ffn=p->XN[IP1];
        result.write((char*)&ffn, sizeof(float));

        ffn=p->YN[JP1];
        result.write((char*)&ffn, sizeof(float));

        ffn=b->bednode(i,j);
        result.write((char*)&ffn, sizeof(float));
    }

    //  Velocities
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol1a(b->P));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->sl_ipol2a(b->Q));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->sl_ipol4(b->ws));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Pressure
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(b->press));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Elevation
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(b->bed));
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
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(b->test));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  Connectivity
    iin=sizeof(int)*(p->polygon_sum)*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
        // Triangle 1
        iin=int(b->nodeval(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(b->nodeval(i,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(b->nodeval(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        // Triangle 2
        iin=int(b->nodeval(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(b->nodeval(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(b->nodeval(i-1,j))-1;
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

    ++printbedcount;
}
