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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>
#include<sys/stat.h>
#include<sys/types.h>

void particle_pls::print(lexer* p, fdm* a, ghostcell* pgc)
{

}

void particle_pls::print_ascii(lexer* p, fdm* a, ghostcell* pgc)
{
//pos
if(posactive-pcount>0)
{
    char name[100];
    sprintf(name,"./REEF3D_PLS/POS-%i-%i.dat",p->count,p->mpirank+1);
	ofstream result;
	result.open(name);

	for(n=0;n<posactive;++n)
	if(posflag[n]>0)
    result<<setprecision(5)<<pos[n][0]+p->originx<<",\t "<<pos[n][1]+p->originy<<",\t "<<pos[n][2]+p->originz<<endl;//",\t "<<pos[n][3]<<",\t "<<pos[n][4]<<endl;//

    result.close();
}

//neg
if(negactive-ncount>0)
{
    char name[100];
    sprintf(name,"./REEF3D_PLS/NEG-%i-%i.dat",p->count,p->mpirank+1);
    ofstream result;
	result.open(name);

    for(n=0;n<negactive;++n)
    if(negflag[n]>0)
    result<<setprecision(5)<<neg[n][0]+p->originx<<",\t "<<neg[n][1]+p->originy<<",\t "<<neg[n][2]+p->originz<<endl;//",\t "<<neg[n][3]<<",\t "<<neg[n][4]<<endl;

    result.close();
}

}

