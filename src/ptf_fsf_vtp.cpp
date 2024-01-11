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

#include"ptf_fsf_vtp.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

ptf_fsf_vtp::ptf_fsf_vtp(lexer* p, fdm_ptf *e, ghostcell *pgc):nodefill(p),vertice(p),nodeflag(p),eta(p),interfac(1.6),zero(0.0)
{
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_FSF",0777);
	
	ptf_fsfprintcount=0;
}

ptf_fsf_vtp::~ptf_fsf_vtp()
{
}

void ptf_fsf_vtp::start(lexer *p, fdm_ptf *e, ghostcell *pgc)
{	
	triangulation(p,e,pgc,e->phi);
	reconstruct(p,e,e->phi);
	
	print(p,e,pgc);
	++ptf_fsfprintcount;
	
	finalize(p,a);
} 

