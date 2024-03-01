/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"sflow_print_wsf.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_vtp_bed::sflow_vtp_bed(lexer *p, fdm2D *b)
{
	if(p->I40==0)
    {
	printbedtime=0.0;
    }
	
	printbedcount=0;
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_SFLOW_VTP_BED",0777);
}

sflow_vtp_bed::~sflow_vtp_bed()
{
}

void sflow_vtp_bed::start(lexer *p, fdm2D* b, ghostcell* pgc, sediment *psed)
{	
	pgc->gcsl_start4(p,b->depth,50);
	pgc->gcsl_start4(p,b->bed,50);
    pgc->gcsl_start4(p,b->test,50);
	
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
	if(p->mpirank==0)
    pvtu(p,b,pgc,psed);
    
	name_iter(p,b,pgc);
    
    b->bed.ggcpol(p);
    
    // bednode upate
    TPSLICELOOP
    {
    pip=4;
    
    b->bednode(i,j) = 0.25*(b->bed(i,j) + b->bed(i+1,j) + b->bed(i,j+1) + b->bed(i,j));

    if(p->flagslice4[Im1Jm1]<0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]>0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j) + b->bed(i,j+1) + b->bed(i,j));
    
    if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]<0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]>0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i,j+1) + b->bed(i,j));
    
    if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]<0 && p->flagslice4[IJ]>0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i+1,j) + b->bed(i,j));
    
    if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]<0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i+1,j) + b->bed(i,j+1));
    
    pip=0;
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
	
	// Open File
	ofstream result;
	result.open(name, ios::binary);
    
    // offsets
    n=0;
	offset[n]=0;
	++n;
	
	// Points
    offset[n]=offset[n-1]+8*(p->pointnum2D)*3+4;
    ++n;
	
	// velocity
	offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
	++n;
	
    // pressure
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // sediment bedlaod
	if(p->P76==1)
	psed->offset_vtp_bedload(p,pgc,result,offset,n);

    // sediment parameters 1
	if(p->P77==1)
	psed->offset_vtp_parameter1(p,pgc,result,offset,n);

    // sediment parameters 2
	if(p->P78==1)
	psed->offset_vtp_parameter2(p,pgc,result,offset,n);

    // bed shear stress
	if(p->P79>=1)
	psed->offset_vtp_bedshear(p,pgc,result,offset,n);
    
    // test
    if(p->P23==1)
    {
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }
	
	// Cells
    offset[n]=offset[n-1] + 4*p->polygon_sum*3+4;
    ++n;
    offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;
	offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;
	
	
	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum2D<<"\" NumberOfPolys=\""<<p->polygon_sum<<"\">"<<endl;
    
    n=0;
	result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float64\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
	
    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
    if(p->P76==1)
	psed->name_vtu_bedload(p,pgc,result,offset,n);
    
    if(p->P77==1)
	psed->name_vtu_parameter1(p,pgc,result,offset,n);

    if(p->P78==1)
	psed->name_vtu_parameter2(p,pgc,result,offset,n);

	if(p->P79>=1)
	psed->name_vtu_bedshear(p,pgc,result,offset,n);
    
    if(p->P23==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    result<<"</PointData>"<<endl;

    

    result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Polys>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</PolyData>"<<endl;
    
    
    //----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";
	
	//  XYZ
	iin=8*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ddn=p->XN[IP1];
	result.write((char*)&ddn, sizeof (double));

	ddn=p->YN[JP1];
	result.write((char*)&ddn, sizeof (double));

	ddn=b->bednode(i,j);
	result.write((char*)&ddn, sizeof (double));
	}
	
    //  Velocities
    iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol1a(b->P));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->sl_ipol2a(b->Q));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->sl_ipol4(b->ws));
	result.write((char*)&ffn, sizeof (float));
	}

	//  Pressure
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->press));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Elevation
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->bed));
	result.write((char*)&ffn, sizeof (float));
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
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->test));
	result.write((char*)&ffn, sizeof (float));
	}
    }

    //  Connectivity
    iin=4*(p->polygon_sum)*3;
    result.write((char*)&iin, sizeof (int));
    SLICEBASELOOP
	{
	// Triangle 1
	iin=int(b->nodeval(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i,j))-1;
	result.write((char*)&iin, sizeof (int));
	
	
	// Triangle 2
	iin=int(b->nodeval(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i,j))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i-1,j))-1;
	result.write((char*)&iin, sizeof (int));
	}
    
    
    //  Offset of Connectivity
    iin=4*(p->polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=(n+1)*3;
	result.write((char*)&iin, sizeof (int));
	}
    
//  Cell types
    iin=4*(p->polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

    result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
	
	++printbedcount;

}


