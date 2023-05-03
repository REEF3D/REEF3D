/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::print_forces_vtp(lexer* p, fdm* a, ghostcell *pgc)
{
	int polygon_num3,polygon_sum3,polygon_sum,vertice_num;
    
	if(p->mpirank==0)
    forces_pvtp(p,a,pgc);
	
	name_iter(p,a,pgc);

	ofstream result;
	result.open(name, ios::binary);
	//---------------------------------------------
	
	polygon_num=facount;
	
	polygon_sum=0;
	for(n=0;n<polygon_num;++n)
	polygon_sum+=numpt[n];
	
	polygon_sum3=polygon_num3=0;
	for(n=0;n<polygon_num;++n)
	if(numpt[n]==4)
	{
	polygon_sum3+=numpt[n];
	++polygon_num3;
	}  
	
	vertice_num = ccptcount;
	
	//---------------------------------------------
    n=0;
	offset[n]=0;
	++n;
    offset[n]=offset[n-1] + 4*(vertice_num)*3 + 4;
    ++n;
	//Data
	offset[n]=offset[n-1] + 4*vertice_num*3+ 4;
    ++n;
	offset[n]=offset[n-1] + 4*vertice_num+ 4;
    ++n;
	//End Data
    offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
    offset[n]=offset[n-1] + 4*polygon_num+ 4;
    ++n;
	offset[n]=offset[n-1] + 4*polygon_num+ 4;
    ++n;
	//---------------------------------------------
	
	

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<vertice_num<<"\" NumberOfPolys=\""<<polygon_num<<"\">"<<endl;

    n=0;
    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</PointData>"<<endl;

    result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	result<<"</Polys>"<<endl;
	

    result<<"</Piece>"<<endl;
    result<<"</PolyData>"<<endl;

//----------------------------------------------------------------------------

    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

//  XYZ
	iin=4*vertice_num*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<vertice_num;++n)
	{
	ffn=ccpt[n][0];
	result.write((char*)&ffn, sizeof (float));

	ffn=ccpt[n][1];
	result.write((char*)&ffn, sizeof (float));

	ffn=ccpt[n][2];
	result.write((char*)&ffn, sizeof (float));
	}
	
//  Velocity
	iin=4*vertice_num*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<vertice_num;++n)
	{
	ffn=float(p->ccipol1(a->u,ccpt[n][0]-p->originx,ccpt[n][1]-p->originy,ccpt[n][2]-p->originz));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ccipol2(a->v,ccpt[n][0]-p->originx,ccpt[n][1]-p->originy,ccpt[n][2]-p->originz));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ccipol3(a->w,ccpt[n][0]-p->originx,ccpt[n][1]-p->originy,ccpt[n][2]-p->originz));
	result.write((char*)&ffn, sizeof (float));
	}
	
	
//  Pressure
	iin=4*vertice_num;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<vertice_num;++n)
	{
            /*i = int((ccpt[n][0]-p->originx)/p->DXM);
            j = int((ccpt[n][1]-p->originy)/p->DXM);
            k = int((ccpt[n][2]-p->originz)/p->DXM);
            
            nx = (a->solid(i+1,j,k)-a->solid(i-1,j,k))/p->DXM;
            ny = (a->solid(i,j+1,k)-a->solid(i,j-1,k))/p->DXM;
            nz = (a->solid(i,j,k+1)-a->solid(i,j,k-1))/p->DXM;
            
            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx/=norm>1.0e-20?norm:1.0e20;
            ny/=norm>1.0e-20?norm:1.0e20;
            nz/=norm>1.0e-20?norm:1.0e20;*/
            
	ffn=float(p->ccipol4(a->press,ccpt[n][0]-p->originx,ccpt[n][1]-p->originy,ccpt[n][2]-p->originz));
	result.write((char*)&ffn, sizeof (float));
	}

//  Connectivity POLYGON
    iin=4*polygon_sum;
    result.write((char*)&iin, sizeof (int));
    for(n=0;n<polygon_num;++n)
	{
		if(numpt[n]==3)
		{
		iin=facet[n][0];
		result.write((char*)&iin, sizeof (int));
		
		iin=facet[n][1];
		result.write((char*)&iin, sizeof (int));
		
		iin=facet[n][2];
		result.write((char*)&iin, sizeof (int));
		}
		
		if(numpt[n]==4)
		{
		iin=facet[n][0];
		result.write((char*)&iin, sizeof (int));
		
		iin=facet[n][1];
		result.write((char*)&iin, sizeof (int));
		
		iin=facet[n][3];
		result.write((char*)&iin, sizeof (int));
		
		iin=facet[n][2];
		result.write((char*)&iin, sizeof (int));
		}
	}

//  Offset of Connectivity
    iin=4*polygon_num;
    result.write((char*)&iin, sizeof (int));
	iin=0;
	for(n=0;n<polygon_num;++n)
	{
	iin+= numpt[n];
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*polygon_num;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<polygon_num;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();	
}


void sixdof_gc::forces_pvtp(lexer* p, fdm* a, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = forceprintcount;

    if(p->P15==2)
    num = p->count;
	
	if(p->P14==0)
	{
    if(num<10)
	sprintf(name,"REEF3D-SOLID-00000%i.pvtp",num);

	if(num<100&&num>9)
	sprintf(name,"REEF3D-SOLID-0000%i.pvtp",num);

	if(num<1000&&num>99)
	sprintf(name,"REEF3D-SOLID-000%i.pvtp",num);

	if(num<10000&&num>999)
	sprintf(name,"REEF3D-SOLID-00%i.pvtp",num);

	if(num<100000&&num>9999)
	sprintf(name,"REEF3D-SOLID-0%i.pvtp",num);

	if(num>99999)
	sprintf(name,"REEF3D-SOLID-%i.pvtp",num);
	}

	if(p->P14==1)
	{
    if(num<10)
	sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00000%i.pvtp",num);

	if(num<100&&num>9)
	sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0000%i.pvtp",num);

	if(num<1000&&num>99)
	sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-000%i.pvtp",num);

	if(num<10000&&num>999)
	sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00%i.pvtp",num);

	if(num<100000&&num>9999)
	sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0%i.pvtp",num);

	if(num>99999)
	sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i.pvtp",num);
	}

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PPolyData  GhostLevel=\"0\">"<<endl;


	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
	
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;
	result<<"</PPointData>"<<endl;
	
	result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"/>"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"/>"<<endl;
	result<<"</Polys>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,a,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PPolyData>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void sixdof_gc::piecename(lexer* p, fdm* a,  ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = forceprintcount;

    if(p->P15==2)
    num = p->count;

	if(n<9)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SOLID-00000%i-0000%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SOLID-0000%i-0000%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SOLID-000%i-0000%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SOLID-00%i-0000%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SOLID-0%i-0000%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SOLID-%i-0000%i.vtp",num,n+1);
	}

	if(n<99&&n>8)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SOLID-00000%i-000%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SOLID-0000%i-000%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SOLID-000%i-000%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SOLID-00%i-000%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SOLID-0%i-000%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SOLID-%i-000%i.vtp",num,n+1);
	}
	if(n<999&&n>98)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SOLID-00000%i-00%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SOLID-0000%i-00%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SOLID-000%i-00%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SOLID-00%i-00%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SOLID-0%i-00%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SOLID-%i-00%i.vtp",num,n+1);
	}

	if(n<9999&&n>998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SOLID-00000%i-0%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SOLID-0000%i-0%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SOLID-000%i-0%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SOLID-00%i-0%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SOLID-0%i-0%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SOLID-%i-0%i.vtp",num,n+1);
	}

	if(n>9998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SOLID-00000%i-%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SOLID-0000%i-%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SOLID-000%i-%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SOLID-00%i-%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SOLID-0%i-%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SOLID-%i-%i.vtp",num,n+1);
	}


}


void sixdof_gc::name_iter(lexer* p,fdm* a,ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = forceprintcount;

    if(p->P15==2)
    num = p->count;

if(p->P14==0)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"REEF3D-SOLID-00000%i-0000%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-SOLID-0000%i-0000%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-SOLID-000%i-0000%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-SOLID-00%i-0000%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-SOLID-0%i-0000%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-SOLID-%i-0000%i.vtp",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"REEF3D-SOLID-00000%i-000%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-SOLID-0000%i-000%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-SOLID-000%i-000%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-SOLID-00%i-000%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-SOLID-0%i-000%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-SOLID-%i-000%i.vtp",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"REEF3D-SOLID-00000%i-00%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-SOLID-0000%i-00%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-SOLID-000%i-00%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-SOLID-00%i-00%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-SOLID-0%i-00%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-SOLID-%i-00%i.vtp",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"REEF3D-SOLID-00000%i-0%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-SOLID-0000%i-0%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-SOLID-000%i-0%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-SOLID-00%i-0%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-SOLID-0%i-0%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-SOLID-%i-0%i.vtp",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"REEF3D-SOLID-00000%i-%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-SOLID-0000%i-%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-SOLID-000%i-%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-SOLID-00%i-%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-SOLID-0%i-%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-SOLID-%i-%i.vtp",num,p->mpirank+1);
	}
}

if(p->P14==1)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00000%i-0000%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0000%i-0000%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-000%i-0000%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00%i-0000%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0%i-0000%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i-0000%i.vtp",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00000%i-000%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0000%i-000%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-000%i-000%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00%i-000%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0%i-000%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i-000%i.vtp",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00000%i-00%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0000%i-00%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-000%i-00%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00%i-00%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0%i-00%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i-00%i.vtp",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00000%i-0%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0000%i-0%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-000%i-0%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00%i-0%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0%i-0%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i-0%i.vtp",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00000%i-%i.vtp",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0000%i-%i.vtp",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-000%i-%i.vtp",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-00%i-%i.vtp",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-0%i-%i.vtp",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i-%i.vtp",num,p->mpirank+1);
	}
}




}




