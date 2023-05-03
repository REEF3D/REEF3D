/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fnpf_vtu3D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fnpf_vtu3D::pvtu(lexer *p, ghostcell* pgc)
{	
	int num=0;
    

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	
	if(p->P14==0)
	{
    if(num<10)
	sprintf(name,"REEF3D-FNPF-00000%i.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"REEF3D-FNPF-0000%i.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"REEF3D-FNPF-000%i.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"REEF3D-FNPF-00%i.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"REEF3D-FNPF-0%i.pvtu",num);

	if(num>99999)
	sprintf(name,"REEF3D-FNPF-%i.pvtu",num);
	}

	if(p->P14==1)
	{
    if(num<10)
	sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-00000%i.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-0000%i.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-000%i.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-00%i.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-0%i.pvtu",num);

	if(num>99999)
	sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-%i.pvtu",num);
	}

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;
	
	
	
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"Fi\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
    if(p->P23==1)
	result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
    if(p->P110==1)
	result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>"<<endl;
    if(p->P25==1)
        result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;
	result<<"</PPointData>"<<endl;
	
    result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
    
	result<<"<Cells>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
	result<<"</Cells>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PUnstructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();

}

void fnpf_vtu3D::piecename(lexer *p, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	if(n<9)
	{
		if(num<10)
		sprintf(pname,"REEF3D-FNPF-00000%i-0000%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-FNPF-0000%i-0000%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-FNPF-000%i-0000%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-FNPF-00%i-0000%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-FNPF-0%i-0000%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-FNPF-%i-0000%i.vtu",num,n+1);
	}

	if(n<99&&n>8)
	{
		if(num<10)
		sprintf(pname,"REEF3D-FNPF-00000%i-000%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-FNPF-0000%i-000%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-FNPF-000%i-000%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-FNPF-00%i-000%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-FNPF-0%i-000%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-FNPF-%i-000%i.vtu",num,n+1);
	}
	if(n<999&&n>98)
	{
		if(num<10)
		sprintf(pname,"REEF3D-FNPF-00000%i-00%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-FNPF-0000%i-00%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-FNPF-000%i-00%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-FNPF-00%i-00%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-FNPF-0%i-00%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-FNPF-%i-00%i.vtu",num,n+1);
	}

	if(n<9999&&n>998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-FNPF-00000%i-0%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-FNPF-0000%i-0%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-FNPF-000%i-0%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-FNPF-00%i-0%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-FNPF-0%i-0%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-FNPF-%i-0%i.vtu",num,n+1);
	}

	if(n>9998)
	{
		if(num<10)
		sprintf(pname,"REEF3D-FNPF-00000%i-%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"REEF3D-FNPF-0000%i-%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"REEF3D-FNPF-000%i-%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"REEF3D-FNPF-00%i-%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"REEF3D-FNPF-0%i-%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"REEF3D-FNPF-%i-%i.vtu",num,n+1);
	}


}
