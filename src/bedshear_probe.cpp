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

#include"bedshear_probe.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

bedshear_probe::bedshear_probe(lexer *p, fdm* a, ghostcell *pgc)
{
    p->Iarray(iloc,p->P125);
	p->Iarray(jloc,p->P125);
	p->Iarray(flag,p->P125);
	p->Darray(bsg,p->P125);
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_SedimentPoint",0777);
	
    if(p->mpirank==0 && p->P125>0)
    {
    // open file
	if(p->P14==0)
    bsgout.open("REEF3D-CFD-Sediment-Bedshear.dat");
	
	if(p->P14==1)
	bsgout.open("./REEF3D_CFD_SedimentPoint/REEF3D-CFD-Sediment-Bedshear.dat");

    bsgout<<"number of gauges:  "<<p->P125<<endl<<endl;
    bsgout<<"x_coord     y_coord"<<endl;
    for(n=0;n<p->P125;++n)
    bsgout<<n+1<<"\t "<<p->P125_x[n]<<"\t "<<p->P125_y[n]<<endl;

    bsgout<<endl<<endl;

    bsgout<<"time";
    for(n=0;n<p->P125;++n)
    bsgout<<"\t P"<<n+1;

    bsgout<<endl<<endl;
    }
	

    ini_location(p,a,pgc);
}

bedshear_probe::~bedshear_probe()
{
    bsgout.close();
}

void bedshear_probe::bedshear_gauge(lexer *p, fdm *a, ghostcell *pgc, sediment *psed)
{
    double zval=0.0;

    for(n=0;n<p->P125;++n)
    bsg[n]=-1.0e20;

	
    for(n=0;n<p->P125;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];

        bsg[n] = psed->bedshear_point(p,a,pgc);
    }
	
	for(n=0;n<p->P121;++n)
    bsg[n]=pgc->globalmax(bsg[n]);

    // write to file
    if(p->mpirank==0)
    {
    bsgout<<p->sedtime<<"\t ";
    for(n=0;n<p->P125;++n)
    bsgout<<bsg[n]<<"  \t  ";
    bsgout<<endl;
    }
}

void bedshear_probe::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void bedshear_probe::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
	
	
    int check;

    for(n=0;n<p->P125;++n)
    {
    iloc[n] = p->posc_i(p->P125_x[n]); 
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n] = p->posc_j(p->P125_y[n]); 

    check=ij_boundcheck(p,a,iloc[n],jloc[n],0);

	
    if(check==1)
    flag[n]=1;
    }
}


int bedshear_probe::conv(double a)
{

int b,c;
double d,diff;

c= int( a);
d=double(c);
diff=a-d;

b=c;

if(diff>0.5)
b=c+1;

if(diff<=-0.5)
b=c-1;

return b;

}
