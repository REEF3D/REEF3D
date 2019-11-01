/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"vtu3D.h"
#include<string>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"heat.h"
#include"vorticity.h"
#include"data.h"
#include"concentration.h"
#include"sediment.h"

void vtu3D::pvtu(fdm* a, lexer* p, ghostcell* pgc, turbulence *pturb, heat *pheat, data *pdata, concentration *pconc, sediment *psed)
{
    int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;
	
	if(p->P14==0)
	{
    if(num<10)
	sprintf(name,"REEF3D-CFD-00000%d.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"REEF3D-CFD-0000%d.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"REEF3D-CFD-000%d.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"REEF3D-CFD-00%d.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"REEF3D-CFD-0%d.pvtu",num);

	if(num>99999)
	sprintf(name,"REEF3D-CFD-%d.pvtu",num);
	}

	if(p->P14==1)
	{
    if(num<10)
	sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00000%d.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0000%d.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-000%d.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00%d.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0%d.pvtu",num);

	if(num>99999)
	sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-%d.pvtu",num);
	}

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;

	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	
	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;

	pturb->name_pvtu(p,a,pgc,result);

	result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>"<<endl;

	pheat->name_vtu(p,a,pgc,result,offset,n);

    pvort->name_vtu(p,a,pgc,result,offset,n);
	
	pdata->name_vtu(p,a,pgc,result,offset,n);
	
	pconc->name_vtu(p,a,pgc,result,offset,n);
    
    if(p->P24==1)
    result<<"<PDataArray type=\"Float32\" Name=\"rho\"/>"<<endl;
    
    if(p->P71==1)
    result<<"<PDataArray type=\"Float32\" Name=\"viscosity\"/>"<<endl;
    
    if(p->P78==1)
    result<<"<PDataArray type=\"Float32\" Name=\"velocity scalar\"/>"<<endl;
    
    if(p->A10==3)
    result<<"<PDataArray type=\"Float32\" Name=\"Fi\"/>"<<endl;
	
	if(p->P26==1)
	{
	result<<"<PDataArray type=\"Float32\" Name=\"cbed\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"conc\"/>"<<endl;
	}
	
	if(p->P27==1)
	result<<"<PDataArray type=\"Float32\" Name=\"topo\"/>"<<endl;
	
	if(p->P79==1)
	psed->name_vtu(p,a,pgc,result,offset,n);
    
    if(p->P23==1)
	result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
	
	result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
	
    if(p->P25==1)
	result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;
    
	if(p->P28==1)
	result<<"<PDataArray type=\"Float32\" Name=\"floating\"/>"<<endl;
	
	if(p->P29==1)
	result<<"<PDataArray type=\"Float32\" Name=\"walldist\"/>"<<endl;
	
	result<<"</PPointData>"<<endl;

	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
	
	result<<"<Cells>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
    ++n;
	result<<"</Cells>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(a,p,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PUnstructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void vtu3D::piecename(fdm* a, lexer* p, ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	if(n<9)
	{
		if(num<10)
		sprintf(pname,"REEF3D-CFD-00000%d-0000%d.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-CFD-0000%d-0000%d.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-CFD-000%d-0000%d.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-CFD-00%d-0000%d.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-CFD-0%d-0000%d.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-CFD-%d-0000%d.vtu",num,n+1);
	}

	if(n<99&&n>8)
	{
		if(num<10)
		sprintf(pname,"REEF3D-CFD-00000%d-000%d.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-CFD-0000%d-000%d.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-CFD-000%d-000%d.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-CFD-00%d-000%d.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-CFD-0%d-000%d.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-CFD-%d-000%d.vtu",num,n+1);
	}
	if(n<999&&n>98)
	{
		if(num<10)
		sprintf(pname,"REEF3D-CFD-00000%d-00%d.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-CFD-0000%d-00%d.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-CFD-000%d-00%d.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-CFD-00%d-00%d.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-CFD-0%d-00%d.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-CFD-%d-00%d.vtu",num,n+1);
	}

	if(n<9999&&n>998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-CFD-00000%d-0%d.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-CFD-0000%d-0%d.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-CFD-000%d-0%d.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-CFD-00%d-0%d.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-CFD-0%d-0%d.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-CFD-%d-0%d.vtu",num,n+1);
	}

	if(n>9998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-CFD-00000%d-%d.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-CFD-0000%d-%d.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-CFD-000%d-%d.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-CFD-00%d-%d.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-CFD-0%d-%d.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-CFD-%d-%d.vtu",num,n+1);
	}


}
