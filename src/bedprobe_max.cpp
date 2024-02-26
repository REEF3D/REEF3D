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

#include"bedprobe_max.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

bedprobe_max::bedprobe_max(lexer *p, fdm* a, ghostcell *pgc)
{
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_SedimentMax",0777);
	
    if(p->mpirank==0 && p->P122>0)
    {
    // open file
	wsfout.open("./REEF3D_CFD_SedimentMax/REEF3D-CFD-Sediment-Max.dat");

    wsfout<<"time  maximum erosion"<<endl<<endl;
    }
}

bedprobe_max::~bedprobe_max()
{
    wsfout.close();
}

void bedprobe_max::bed_max(lexer *p, fdm *a, ghostcell *pgc)
{
    // write to file
    if(p->mpirank==0)
    wsfout<<p->sedtime<<" \t "<<p->bedmin<<endl;
}
