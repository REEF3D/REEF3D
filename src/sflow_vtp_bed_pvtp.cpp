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
#include"sediment.h"

void sflow_vtp_bed::pvtu(lexer *p, fdm2D* b, ghostcell* pgc, sediment *psed)
{	
	int num=0;

    if(p->P15==1)
    num = printbedcount;

    if(p->P15==2)
    num = p->count;
	
	if(p->P14==0)
	{
    if(num<10)
	sprintf(name,"REEF3D-SFLOW-BED-00000%i.pvtp",num);

	if(num<100&&num>9)
	sprintf(name,"REEF3D-SFLOW-BED-0000%i.pvtp",num);

	if(num<1000&&num>99)
	sprintf(name,"REEF3D-SFLOW-BED-000%i.pvtp",num);

	if(num<10000&&num>999)
	sprintf(name,"REEF3D-SFLOW-BED-00%i.pvtp",num);

	if(num<100000&&num>9999)
	sprintf(name,"REEF3D-SFLOW-BED-0%i.pvtp",num);

	if(num>99999)
	sprintf(name,"REEF3D-SFLOW-BED-%i.pvtp",num);
	}

	if(p->P14==1)
	{
    if(num<10)
	sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-00000%i.pvtp",num);

	if(num<100&&num>9)
	sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-0000%i.pvtp",num);

	if(num<1000&&num>99)
	sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-000%i.pvtp",num);

	if(num<10000&&num>999)
	sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-00%i.pvtp",num);

	if(num<100000&&num>9999)
	sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-0%i.pvtp",num);

	if(num>99999)
	sprintf(name,"./REEF3D_SFLOW_VTP_BED/REEF3D-SFLOW-BED-%i.pvtp",num);
	}

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PPolyData  GhostLevel=\"0\">"<<endl;
	
	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
	
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
    
    if(p->P76==1)
	psed->name_pvtu_bedload(p,pgc,result);
    
    if(p->P77==1)
	psed->name_pvtu_parameter1(p,pgc,result);

    if(p->P78==1)
	psed->name_pvtu_parameter2(p,pgc,result);

	if(p->P79>=1)
	psed->name_pvtu_bedshear(p,pgc,result);
    
    if(p->P23==1)
    result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
	result<<"</PPointData>"<<endl;
	
	result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
	result<<"</Polys>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,b,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PPolyData>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();

}

void sflow_vtp_bed::piecename(lexer *p, fdm2D *b, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printbedcount;

    if(p->P15==2)
    num = p->count;

	if(n<9)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SFLOW-BED-00000%i-0000%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SFLOW-BED-0000%i-0000%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SFLOW-BED-000%i-0000%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SFLOW-BED-00%i-0000%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SFLOW-BED-0%i-0000%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SFLOW-BED-%i-0000%i.vtp",num,n+1);
	}

	if(n<99&&n>8)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SFLOW-BED-00000%i-000%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SFLOW-BED-0000%i-000%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SFLOW-BED-000%i-000%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SFLOW-BED-00%i-000%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SFLOW-BED-0%i-000%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SFLOW-BED-%i-000%i.vtp",num,n+1);
	}
	if(n<999&&n>98)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SFLOW-BED-00000%i-00%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SFLOW-BED-0000%i-00%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SFLOW-BED-000%i-00%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SFLOW-BED-00%i-00%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SFLOW-BED-0%i-00%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SFLOW-BED-%i-00%i.vtp",num,n+1);
	}

	if(n<9999&&n>998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SFLOW-BED-00000%i-0%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SFLOW-BED-0000%i-0%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SFLOW-BED-000%i-0%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SFLOW-BED-00%i-0%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SFLOW-BED-0%i-0%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SFLOW-BED-%i-0%i.vtp",num,n+1);
	}

	if(n>9998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-SFLOW-BED-00000%i-%i.vtp",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-SFLOW-BED-0000%i-%i.vtp",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-SFLOW-BED-000%i-%i.vtp",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-SFLOW-BED-00%i-%i.vtp",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-SFLOW-BED-0%i-%i.vtp",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-SFLOW-BED-%i-%i.vtp",num,n+1);
	}


}
