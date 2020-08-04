/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_vtp_bed.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_print_wsf.h"
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
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SFLOW_VTP_BED",0777);
	
	
}

sflow_vtp_bed::~sflow_vtp_bed()
{
}

void sflow_vtp_bed::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	pgc->gcsl_start4(p,b->depth,50);
	pgc->gcsl_start4(p,b->bed,50);
    pgc->gcsl_start4(p,b->test,50);
	
	// Print out based on iteration
    if((p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)  || (p->count==0 &&  p->P30<0.0))
    {
    print2D(p,b,pgc);
    }
		
    // Print out based on time
    if((p->simtime>printbedtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
    {
    print2D(p,b,pgc);
		
    printbedtime+=p->P30;
    }
}

void sflow_vtp_bed::print2D(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	if(p->mpirank==0)
    pvtu(p,b,pgc);
    
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
    offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
    ++n;
	
	// velocity
	offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
	++n;
	
	// bedload
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // bedchange
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	
    // pressure
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	
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
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
	
    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"waterlevel\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"bedload\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"bedchange\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
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
	iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
    
	ffn=float(float(i+1)*p->DXM+p->originx);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(float(j+1)*p->DXM+p->originy);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(b->bednode(i,j));
	result.write((char*)&ffn, sizeof (float));
	}
	
    //  Velocities
    iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(pgc->gcsl_ipol1a(p,b->hx));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(pgc->gcsl_ipol2a(p,b->hy));
	result.write((char*)&ffn, sizeof (float));
	
	ffn=float(pgc->gcsl_ipol4(p,b->hp));
	result.write((char*)&ffn, sizeof (float));
	}

    //  bedload
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(pgc->gcsl_ipol4(p,b->qb));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  bedchange
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(pgc->gcsl_ipol4(p,b->zb));
	result.write((char*)&ffn, sizeof (float));
	}
	
	//  Test
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(pgc->gcsl_ipol4(p,b->test));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Elevation
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(pgc->gcsl_ipol4(p,b->bed));
	result.write((char*)&ffn, sizeof (float));
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


