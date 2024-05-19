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

#include"sflow_vtp_fsf.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<string>

void sflow_vtp_fsf::name_iter(lexer *p, fdm2D* b, ghostcell* pgc)
{	
    int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;


    sprintf(name,"./REEF3D_SFLOW_VTP_FSF/REEF3D-SFLOW-FSF-%08i-%06i.vtp",num,p->mpirank+1);
}

