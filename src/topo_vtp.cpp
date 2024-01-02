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

#include"topo_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

topo_vtp::topo_vtp(lexer* p, fdm *a, ghostcell *pgc)
{
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_Topo",0777);
	
	topoprintcount=0;
    polygon_sum=0;
    
    SLICEBASELOOP
	++polygon_sum;
    
    polygon_sum*=2;
}

topo_vtp::~topo_vtp()
{
}

void topo_vtp::start(lexer *p, fdm *a, ghostcell *pgc, sediment *psed)
{	
	print(p,a,pgc,psed);
	++topoprintcount;
} 

