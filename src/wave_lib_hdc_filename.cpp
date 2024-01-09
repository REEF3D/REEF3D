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

#include"wave_lib_hdc.h"
#include"lexer.h"

void wave_lib_hdc::filename_single(lexer *p, ghostcell *pgc,int num)
{
    if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000000%i-0000%i.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000000%i-0000%i.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00000%i-0000%i.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000%i-0000%i.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000%i-0000%i.r3d",num,p->mpirank+1);

		if(num<1000000&&num>99999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00%i-0000%i.r3d",num,p->mpirank+1);
        
        if(num<10000000&&num>999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0%i-0000%i.r3d",num,p->mpirank+1);
        
        if(num>9999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%i-0000%i.r3d",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000000%i-000%i.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000000%i-000%i.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00000%i-000%i.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000%i-000%i.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000%i-000%i.r3d",num,p->mpirank+1);

		if(num<1000000&&num>99999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00%i-000%i.r3d",num,p->mpirank+1);
        
        if(num<10000000&&num>999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0%i-000%i.r3d",num,p->mpirank+1);
        
        if(num>9999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%i-000%i.r3d",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000000%i-00%i.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000000%i-00%i.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00000%i-00%i.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000%i-00%i.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000%i-00%i.r3d",num,p->mpirank+1);

		if(num<1000000&&num>99999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00%i-00%i.r3d",num,p->mpirank+1);
        
        if(num<10000000&&num>999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0%i-00%i.r3d",num,p->mpirank+1);
        
        if(num>9999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%i-00%i.r3d",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000000%i-0%i.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000000%i-0%i.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00000%i-0%i.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000%i-0%i.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000%i-0%i.r3d",num,p->mpirank+1);

		if(num<1000000&&num>99999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00%i-0%i.r3d",num,p->mpirank+1);
        
        if(num<10000000&&num>999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0%i-0%i.r3d",num,p->mpirank+1);
        
        if(num>9999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%i-0%i.r3d",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000000%i-%i.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000000%i-%i.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00000%i-%i.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000%i-%i.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000%i-%i.r3d",num,p->mpirank+1);

		if(num<1000000&&num>99999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00%i-%i.r3d",num,p->mpirank+1);
        
        if(num<10000000&&num>999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0%i-%i.r3d",num,p->mpirank+1);
        
        if(num>9999999)
		sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%i-%i.r3d",num,p->mpirank+1);
	}
    
    
}

void wave_lib_hdc::filename_continuous(lexer *p, ghostcell *pgc)
{
    if(p->mpirank<9)
    sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0000%i.r3d",p->mpirank+1);

	if(p->mpirank<99&&p->mpirank>8)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-000%i.r3d",p->mpirank+1);

	if(p->mpirank<999&&p->mpirank>98)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-00%i.r3d",p->mpirank+1);

	if(p->mpirank<9999&&p->mpirank>998)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-0%i.r3d",p->mpirank+1);

	if(p->mpirank>9998)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%i.r3d",p->mpirank+1);
 
}

void wave_lib_hdc::filename_header(lexer *p, ghostcell *pgc)
{

    if(p->mpirank<9)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-Header-0000%i.r3d",p->mpirank+1);

	if(p->mpirank<99&&p->mpirank>8)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-Header-000%i.r3d",p->mpirank+1);
    
	if(p->mpirank<999&&p->mpirank>98)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-Header-00%i.r3d",p->mpirank+1);

	if(p->mpirank<9999&&p->mpirank>998)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-Header-0%i.r3d",p->mpirank+1);

	if(p->mpirank>9998)
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-Header-%i.r3d",p->mpirank+1);   
	
}
