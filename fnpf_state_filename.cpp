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

#include"fnpf_state.h"
#include"lexer.h"

void fnpf_state::filename(lexer *p, fdm_fnpf *c, ghostcell *pgc, int num)
{
    
if(p->P14==0)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"REEF3D_FNPF-State-00000%d-0000%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D_FNPF-State-0000%d-0000%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D_FNPF-State-000%d-0000%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D_FNPF-State-00%d-0000%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D_FNPF-State-0%d-0000%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D_FNPF-State-%d-0000%d.r3d",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"REEF3D_FNPF-State-00000%d-000%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D_FNPF-State-0000%d-000%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D_FNPF-State-000%d-000%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D_FNPF-State-00%d-000%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D_FNPF-State-0%d-000%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D_FNPF-State-%d-000%d.r3d",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"REEF3D_FNPF-State-00000%d-00%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D_FNPF-State-0000%d-00%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D_FNPF-State-000%d-00%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D_FNPF-State-00%d-00%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D_FNPF-State-0%d-00%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D_FNPF-State-%d-00%d.r3d",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"REEF3D_FNPF-State-00000%d-0%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D_FNPF-State-0000%d-0%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D_FNPF-State-000%d-0%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D_FNPF-State-00%d-0%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D_FNPF-State-0%d-0%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D_FNPF-State-%d-0%d.r3d",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"REEF3D_FNPF-State-00000%d-%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"REEF3D_FNPF-State-0000%d-%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D_FNPF-State-000%d-%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D_FNPF-State-00%d-%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D_FNPF-State-0%d-%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"REEF3D_FNPF-State-%d-%d.r3d",num,p->mpirank+1);
	}
}

if(p->P14==1)
{
	if(p->mpirank<9)
	{
		if(num<10)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00000%d-0000%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0000%d-0000%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-000%d-0000%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00%d-0000%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0%d-0000%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-%d-0000%d.r3d",num,p->mpirank+1);
	}

	if(p->mpirank<99&&p->mpirank>8)
	{
		if(num<10)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00000%d-000%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0000%d-000%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-000%d-000%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00%d-000%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0%d-000%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-%d-000%d.r3d",num,p->mpirank+1);
	}
	if(p->mpirank<999&&p->mpirank>98)
	{
		if(num<10)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00000%d-00%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0000%d-00%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-000%d-00%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00%d-00%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0%d-00%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-%d-00%d.r3d",num,p->mpirank+1);
	}

	if(p->mpirank<9999&&p->mpirank>998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00000%d-0%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0000%d-0%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-000%d-0%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00%d-0%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0%d-0%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-%d-0%d.r3d",num,p->mpirank+1);
	}

	if(p->mpirank>9998)
	{
		if(num<10)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00000%d-%d.r3d",num,p->mpirank+1);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0000%d-%d.r3d",num,p->mpirank+1);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-000%d-%d.r3d",num,p->mpirank+1);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-00%d-%d.r3d",num,p->mpirank+1);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-0%d-%d.r3d",num,p->mpirank+1);

		if(num>99999)
		sprintf(name,"./REEF3D_FNPF_State/REEF3D_FNPF-State-%d-%d.r3d",num,p->mpirank+1);
	}
}
}

void fnpf_state::filename_header(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->P14==0)
    {
	if(p->mpirank<9)
	sprintf(name,"REEF3D-FNPF-State-Header-0000%d.r3d",p->mpirank+1);

	if(p->mpirank<99&&p->mpirank>8)
	sprintf(name,"REEF3D-FNPF-State-Header-000%d.r3d",p->mpirank+1);
    
	if(p->mpirank<999&&p->mpirank>98)
	sprintf(name,"REEF3D-FNPF-State-Header-00%d.r3d",p->mpirank+1);

	if(p->mpirank<9999&&p->mpirank>998)
	sprintf(name,"REEF3D-FNPF-State-Header-0%d.r3d",p->mpirank+1);

	if(p->mpirank>9998)
	sprintf(name,"REEF3D-FNPF-State-Header-%d.r3d",p->mpirank+1);
    }
    
    if(p->P14==1)
    {
	if(p->mpirank<9)
	sprintf(name,"./REEF3D_FNPF_State/REEF3D-FNPF-State-Header-0000%d.r3d",p->mpirank+1);

	if(p->mpirank<99&&p->mpirank>8)
	sprintf(name,"./REEF3D_FNPF_State/REEF3D-FNPF-State-Header-000%d.r3d",p->mpirank+1);
    
	if(p->mpirank<999&&p->mpirank>98)
	sprintf(name,"./REEF3D_FNPF_State/REEF3D-FNPF-State-Header-00%d.r3d",p->mpirank+1);

	if(p->mpirank<9999&&p->mpirank>998)
	sprintf(name,"./REEF3D_FNPF_State/REEF3D-FNPF-State-Header-0%d.r3d",p->mpirank+1);

	if(p->mpirank>9998)
	sprintf(name,"./REEF3D_FNPF_State/REEF3D-FNPF-State-Header-%d.r3d",p->mpirank+1);
    }
}



