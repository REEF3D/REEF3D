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

#include"bedshear_max.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

bedshear_max::bedshear_max(lexer *p, fdm* a, ghostcell *pgc)
{
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_SedimentMax",0777);
	
    if(p->mpirank==0 && p->P126>0)
    {
    // open file
	if(p->P14==0)
    bsgout.open("REEF3D-CFD-Sediment-Bedshear-Max.dat");
	
	if(p->P14==1)
	bsgout.open("./REEF3D_CFD_SedimentPoint/REEF3D-CFD-Sediment-Bedshear-Max.dat");


    bsgout<<"time";
    bsgout<<"\t  bedshear max";

    bsgout<<endl<<endl;
    }
	

}

bedshear_max::~bedshear_max()
{
    bsgout.close();
}

void bedshear_max::bedshear_maxval(lexer *p, fdm *a, ghostcell *pgc, sediment *psed)
{
    double maxval;

    maxval=-1.0e20;

	
    ILOOP
    JLOOP
    maxval = MAX(maxval, psed->bedshear_point(p,a,pgc));

	
    maxval=pgc->globalmax(maxval);

    // write to file
    if(p->mpirank==0)
    {
    bsgout<<p->sedtime<<"\t ";
    bsgout<<maxval<<endl;
    }
}


