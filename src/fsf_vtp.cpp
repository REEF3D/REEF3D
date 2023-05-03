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

#include"fsf_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

fsf_vtp::fsf_vtp(lexer* p, fdm *a, ghostcell *pgc):nodefill(p),vertice(p),nodeflag(p),eta(p),interfac(1.6),zero(0.0)
{
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_FSF",0777);
	
	fsfprintcount=0;
}

fsf_vtp::~fsf_vtp()
{
}

void fsf_vtp::start(lexer *p, fdm *a, ghostcell *pgc)
{	
	triangulation(p,a,pgc,a->phi);
	reconstruct(p,a,a->phi);
	
	print(p,a,pgc);
	++fsfprintcount;
	
	finalize(p,a);
} 

