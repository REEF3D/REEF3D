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
#include"fnpf_vtu3D.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<string>

void fnpf_vtu3D::name_iter(lexer *p, ghostcell* pgc)
{	
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

    sprintf(name,"./REEF3D_FNPF_VTU/REEF3D-FNPF-%08i-%06i.vtu",num,p->mpirank+1);

}

