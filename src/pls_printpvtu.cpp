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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void particle_pls::pvtu_pos(fdm* a, lexer* p, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	if(p->P14==0)
	{
    if(num<10)
	sprintf(name,"XPLS-POS-00000%i.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"XPLS-POS-0000%i.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"XPLS-POS-000%i.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"XPLS-POS-00%i.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"XPLS-POS-0%i.pvtu",num);

	if(num>99999)
	sprintf(name,"XPLS-POS-%i.pvtu",num);
	}
	
	if(p->P14==1)
	{
	if(num<10)
	sprintf(name,"./REEF3D_PLS/XPLS-POS-00000%i.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"./REEF3D_PLS/XPLS-POS-0000%i.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"./REEF3D_PLS/XPLS-POS-000%i.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"./REEF3D_PLS/XPLS-POS-00%i.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"./REEF3D_PLS/XPLS-POS-0%i.pvtu",num);

	if(num>99999)
	sprintf(name,"./REEF3D_PLS/XPLS-POS-%i.pvtu",num);
	}

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;

	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"radius\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"correction\"/>"<<endl;
	result<<"</PPointData>"<<endl;

	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename_pos(a,p,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PUnstructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}


void particle_pls::pvtu_neg(fdm* a, lexer* p, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	
	if(p->P14==0)
	{
    if(num<10)
	sprintf(name,"XPLS-NEG-00000%i.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"XPLS-NEG-0000%i.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"XPLS-NEG-000%i.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"XPLS-NEG-00%i.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"XPLS-NEG-0%i.pvtu",num);

	if(num>99999)
	sprintf(name,"XPLS-NEG-%i.pvtu",num);
	}
	
	if(p->P14==1)
	{
    if(num<10)
	sprintf(name,"./REEF3D_PLS/XPLS-NEG-00000%i.pvtu",num);

	if(num<100&&num>9)
	sprintf(name,"./REEF3D_PLS/XPLS-NEG-0000%i.pvtu",num);

	if(num<1000&&num>99)
	sprintf(name,"./REEF3D_PLS/XPLS-NEG-000%i.pvtu",num);

	if(num<10000&&num>999)
	sprintf(name,"./REEF3D_PLS/XPLS-NEG-00%i.pvtu",num);

	if(num<100000&&num>9999)
	sprintf(name,"./REEF3D_PLS/XPLS-NEG-0%i.pvtu",num);

	if(num>99999)
	sprintf(name,"./REEF3D_PLS/XPLS-NEG-%i.pvtu",num);
	}
	
	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;

	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"radius\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"correction\"/>"<<endl;
	result<<"</PPointData>"<<endl;

	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename_neg(a,p,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PUnstructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void particle_pls::piecename_pos(fdm* a, lexer* p, ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;


	if(n<9)
	{
		if(num<10)
		sprintf(pname,"XPLS-POS-00000%i-0000%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"XPLS-POS-0000%i-0000%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"XPLS-POS-000%i-0000%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"XPLS-POS-00%i-0000%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"XPLS-POS-0%i-0000%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"XPLS-POS-%i-0000%i.vtu",num,n+1);
	}

	if(n<99&&n>8)
	{
		if(num<10)
		sprintf(pname,"XPLS-POS-00000%i-000%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"XPLS-POS-0000%i-000%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"XPLS-POS-000%i-000%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"XPLS-POS-00%i-000%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"XPLS-POS-0%i-000%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"XPLS-POS-%i-000%i.vtu",num,n+1);
	}
	if(n<999&&n>98)
	{
		if(num<10)
		sprintf(pname,"XPLS-POS-00000%i-00%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"XPLS-POS-0000%i-00%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"XPLS-POS-000%i-00%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"XPLS-POS-00%i-00%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"XPLS-POS-0%i-00%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"XPLS-POS-%i-00%i.vtu",num,n+1);
	}

	if(n<9999&&n>998)
	{
		if(num<10)
		sprintf(pname,"XPLS-POS-00000%i-0%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"XPLS-POS-0000%i-0%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"XPLS-POS-000%i-0%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"XPLS-POS-00%i-0%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"XPLS-POS-0%i-0%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"XPLS-POS-%i-0%i.vtu",num,n+1);
	}

	if(n>9998)
	{
		if(num<10)
		sprintf(pname,"XPLS-POS-00000%i-%i.vtu",num,n+1);

		if(num<100&&num>9)
		sprintf(pname,"XPLS-POS-0000%i-%i.vtu",num,n+1);

		if(num<1000&&num>99)
		sprintf(pname,"XPLS-POS-000%i-%i.vtu",num,n+1);

		if(num<10000&&num>999)
		sprintf(pname,"XPLS-POS-00%i-%i.vtu",num,n+1);

		if(num<100000&&num>9999)
		sprintf(pname,"XPLS-POS-0%i-%i.vtu",num,n+1);

		if(num>99999)
		sprintf(pname,"XPLS-POS-%i-%i.vtu",num,n+1);
	}

}



void particle_pls::piecename_neg(fdm* a, lexer* p, ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

if(n<9)
{
	if(num<10)
	sprintf(pname,"XPLS-NEG-00000%i-0000%i.vtu",num,n+1);

	if(num<100&&num>9)
	sprintf(pname,"XPLS-NEG-0000%i-0000%i.vtu",num,n+1);

	if(num<1000&&num>99)
	sprintf(pname,"XPLS-NEG-000%i-0000%i.vtu",num,n+1);

	if(num<10000&&num>999)
	sprintf(pname,"XPLS-NEG-00%i-0000%i.vtu",num,n+1);

	if(num<100000&&num>9999)
	sprintf(pname,"XPLS-NEG-0%i-0000%i.vtu",num,n+1);

	if(num>99999)
	sprintf(pname,"XPLS-NEG-%i-0000%i.vtu",num,n+1);
}

if(n<99&&n>8)
{
	if(num<10)
	sprintf(pname,"XPLS-NEG-00000%i-000%i.vtu",num,n+1);

	if(num<100&&num>9)
	sprintf(pname,"XPLS-NEG-0000%i-000%i.vtu",num,n+1);

	if(num<1000&&num>99)
	sprintf(pname,"XPLS-NEG-000%i-000%i.vtu",num,n+1);

	if(num<10000&&num>999)
	sprintf(pname,"XPLS-NEG-00%i-000%i.vtu",num,n+1);

	if(num<100000&&num>9999)
	sprintf(pname,"XPLS-NEG-0%i-000%i.vtu",num,n+1);

	if(num>99999)
	sprintf(pname,"XPLS-NEG-%i-000%i.vtu",num,n+1);
}
if(n<999&&n>98)
{
	if(num<10)
	sprintf(pname,"XPLS-NEG-00000%i-00%i.vtu",num,n+1);

	if(num<100&&num>9)
	sprintf(pname,"XPLS-NEG-0000%i-00%i.vtu",num,n+1);

	if(num<1000&&num>99)
	sprintf(pname,"XPLS-NEG-000%i-00%i.vtu",num,n+1);

	if(num<10000&&num>999)
	sprintf(pname,"XPLS-NEG-00%i-00%i.vtu",num,n+1);

	if(num<100000&&num>9999)
	sprintf(pname,"XPLS-NEG-0%i-00%i.vtu",num,n+1);

	if(num>99999)
	sprintf(pname,"XPLS-NEG-%i-00%i.vtu",num,n+1);
}

if(n<9999&&n>998)
{
	if(num<10)
	sprintf(pname,"XPLS-NEG-00000%i-0%i.vtu",num,n+1);

	if(num<100&&num>9)
	sprintf(pname,"XPLS-NEG-0000%i-0%i.vtu",num,n+1);

	if(num<1000&&num>99)
	sprintf(pname,"XPLS-NEG-000%i-0%i.vtu",num,n+1);

	if(num<10000&&num>999)
	sprintf(pname,"XPLS-NEG-00%i-0%i.vtu",num,n+1);

	if(num<100000&&num>9999)
	sprintf(pname,"XPLS-NEG-0%i-0%i.vtu",num,n+1);

	if(num>99999)
	sprintf(pname,"XPLS-NEG-%i-0%i.vtu",num,n+1);
}

if(n>9998)
{
	if(num<10)
	sprintf(pname,"XPLS-NEG-00000%i-%i.vtu",num,n+1);

	if(num<100&&num>9)
	sprintf(pname,"XPLS-NEG-0000%i-%i.vtu",num,n+1);

	if(num<1000&&num>99)
	sprintf(pname,"XPLS-NEG-000%i-%i.vtu",num,n+1);

	if(num<10000&&num>999)
	sprintf(pname,"XPLS-NEG-00%i-%i.vtu",num,n+1);

	if(num<100000&&num>9999)
	sprintf(pname,"XPLS-NEG-0%i-%i.vtu",num,n+1);

	if(num>99999)
	sprintf(pname,"XPLS-NEG-%i-%i.vtu",num,n+1);
}

}


void particle_pls::header_pos(fdm* a,lexer* p,ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

if(p->P14==0)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"XPLS-POS-00000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-POS-0000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-POS-000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-POS-00%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-POS-0%i-0000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-POS-%i-0000%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"XPLS-POS-00000%i-000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-POS-0000%i-000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-POS-000%i-000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-POS-00%i-000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-POS-0%i-000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-POS-%i-000%i.vtu",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"XPLS-POS-00000%i-00%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-POS-0000%i-00%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-POS-000%i-00%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-POS-00%i-00%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-POS-0%i-00%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-POS-%i-00%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"XPLS-POS-00000%i-0%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-POS-0000%i-0%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-POS-000%i-0%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-POS-00%i-0%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-POS-0%i-0%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-POS-%i-0%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"XPLS-POS-00000%i-%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-POS-0000%i-%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-POS-000%i-%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-POS-00%i-%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-POS-0%i-%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-POS-%i-%i.vtu",num,p->mpirank+1);
	}
}


if(p->P14==1)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0%i-0000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-%i-0000%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00000%i-000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0000%i-000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-000%i-000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00%i-000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0%i-000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-%i-000%i.vtu",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00000%i-00%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0000%i-00%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-000%i-00%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00%i-00%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0%i-00%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-%i-00%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00000%i-0%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0000%i-0%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-000%i-0%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00%i-0%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0%i-0%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-%i-0%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00000%i-%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0000%i-%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-000%i-%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-00%i-%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-0%i-%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-POS-%i-%i.vtu",num,p->mpirank+1);
	}
}

}


void particle_pls::header_neg(fdm* a,lexer* p,ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

if(p->P14==0)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"XPLS-NEG-00000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-NEG-0000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-NEG-000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-NEG-00%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-NEG-0%i-0000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-NEG-%i-0000%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"XPLS-NEG-00000%i-000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-NEG-0000%i-000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-NEG-000%i-000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-NEG-00%i-000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-NEG-0%i-000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-NEG-%i-000%i.vtu",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"XPLS-NEG-00000%i-00%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-NEG-0000%i-00%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-NEG-000%i-00%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-NEG-00%i-00%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-NEG-0%i-00%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-NEG-%i-00%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"XPLS-NEG-00000%i-0%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-NEG-0000%i-0%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-NEG-000%i-0%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-NEG-00%i-0%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-NEG-0%i-0%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-NEG-%i-0%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"XPLS-NEG-00000%i-%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"XPLS-NEG-0000%i-%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"XPLS-NEG-000%i-%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"XPLS-NEG-00%i-%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"XPLS-NEG-0%i-%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"XPLS-NEG-%i-%i.vtu",num,p->mpirank+1);
	}
}

if(p->P14==1)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0%i-0000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-%i-0000%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00000%i-000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0000%i-000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-000%i-000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00%i-000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0%i-000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-%i-000%i.vtu",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00000%i-00%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0000%i-00%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-000%i-00%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00%i-00%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0%i-00%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-%i-00%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00000%i-0%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0000%i-0%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-000%i-0%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00%i-0%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0%i-0%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-%i-0%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00000%i-%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0000%i-%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-000%i-%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-00%i-%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-0%i-%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_PLS/XPLS-NEG-%i-%i.vtu",num,p->mpirank+1);
	}
}

}
