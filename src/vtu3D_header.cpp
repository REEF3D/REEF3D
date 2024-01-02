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

#include"vtu3D.h"
#include<string>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vtu3D::name_iter(fdm* a,lexer* p,ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

if(p->P14==0)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"REEF3D-CFD-00000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-CFD-0000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-CFD-000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-CFD-00%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-CFD-0%i-0000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-CFD-%i-0000%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"REEF3D-CFD-00000%i-000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-CFD-0000%i-000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-CFD-000%i-000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-CFD-00%i-000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-CFD-0%i-000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-CFD-%i-000%i.vtu",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"REEF3D-CFD-00000%i-00%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-CFD-0000%i-00%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-CFD-000%i-00%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-CFD-00%i-00%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-CFD-0%i-00%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-CFD-%i-00%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"REEF3D-CFD-00000%i-0%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-CFD-0000%i-0%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-CFD-000%i-0%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-CFD-00%i-0%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-CFD-0%i-0%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-CFD-%i-0%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"REEF3D-CFD-00000%i-%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-CFD-0000%i-%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-CFD-000%i-%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-CFD-00%i-%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-CFD-0%i-%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D-CFD-%i-%i.vtu",num,p->mpirank+1);
	}
}

if(p->P14==1)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-000%i-0000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00%i-0000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0%i-0000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-%i-0000%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00000%i-000%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0000%i-000%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-000%i-000%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00%i-000%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0%i-000%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-%i-000%i.vtu",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00000%i-00%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0000%i-00%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-000%i-00%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00%i-00%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0%i-00%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-%i-00%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00000%i-0%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0000%i-0%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-000%i-0%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00%i-0%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0%i-0%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-%i-0%i.vtu",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00000%i-%i.vtu",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0000%i-%i.vtu",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-000%i-%i.vtu",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-00%i-%i.vtu",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-0%i-%i.vtu",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_CFD_VTU/REEF3D-CFD-%i-%i.vtu",num,p->mpirank+1);
	}
}

    sprintf(epsvar,"epsilon");

	if(p->T10==2||p->T10==12 || p->T10==22)
	sprintf(epsvar,"omega");


}

void vtu3D::name_time(fdm* a,lexer* p,ghostcell* pgc)
{

}

void vtu3D::header(fdm* a,lexer* p,ghostcell* pgc)
{

}
