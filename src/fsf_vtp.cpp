/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"plic_distance.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

fsf_vtp::fsf_vtp(lexer* p, fdm *a, ghostcell *pgc):vertice(p),nodeflag(p),eta(p),interfac(1.6),zero(0.0)
{
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_FSF",0777);
	
	fsfprintcount=0;

    if(p->F80==4)
    pplicd = std::make_shared<plic_distance>(p);
}

fsf_vtp::~fsf_vtp()
{
}

void fsf_vtp::start(lexer *p, fdm *a, ghostcell *pgc)
{	
	// VOF-PLIC: phi is a phase marker; the surface is cut from the distance to the PLIC planes
    field &ls = pplicd ? static_cast<field&>(pplicd->start(p,pgc,a->vof)) : static_cast<field&>(a->phi);

	triangulation(p,a,pgc,ls);
	reconstruct(p,a,ls);
	
	print(p,a);
	++fsfprintcount;
	
	finalize(p,a);
} 

