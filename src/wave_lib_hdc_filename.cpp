/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
    sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%08i-%06i.r3d",num,p->mpirank+1);    
}

void wave_lib_hdc::filename_continuous(lexer *p, ghostcell *pgc)
{
    sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-%06i.r3d",p->mpirank+1);
}

void wave_lib_hdc::filename_header(lexer *p, ghostcell *pgc)
{
	sprintf(name,"./REEF3D_CFD_HDC_Input/REEF3D-HDC-Input-Header-%06i.r3d",p->mpirank+1);  
}
