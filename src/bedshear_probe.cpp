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

#include"bedshear_probe.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<stdio.h>

bedshear_probe::bedshear_probe(lexer *p, ghostcell *pgc)
{
    p->Iarray(iloc,p->P125);
	p->Iarray(jloc,p->P125);
	p->Iarray(flag,p->P125);
	p->Darray(bsg,p->P125);
	
	// Create Folder
	if(p->mpirank==0)
    {
        char folder[40];
        if(p->A10==5)
            snprintf(folder,sizeof(folder),"./REEF3D_NHFLOW_Sediment");
        else
            snprintf(folder,sizeof(folder),"./REEF3D_CFD_Sediment");
	    mkdir(folder,0777);
    }
	
    if(p->mpirank==0 && p->P125>0)
    {
    // open file
        char file[100];
        if(p->A10==5)
            snprintf(file,sizeof(file),"./REEF3D_NHFLOW_Sediment/REEF3D-NHFLOW-Sediment-Bedshear.dat");
        else
            snprintf(file,sizeof(file),"./REEF3D_CFD_Sediment/REEF3D-CFD-Sediment-Bedshear.dat");
	    bsgout.open(file);

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
	

    ini_location(p,pgc);
}

bedshear_probe::~bedshear_probe()
{
    bsgout.close();
}

void bedshear_probe::bedshear_gauge(lexer *p, ghostcell *pgc, sediment *psed)
{
    for(n=0;n<p->P125;++n)
        bsg[n]=-1.0e20;

	
    for(n=0;n<p->P125;++n)
        if(flag[n]>0)
        {

            i=iloc[n];
            j=jloc[n];
            bsg[n] = psed->bedshear_point(p,pgc);
        }
	
	for(n=0;n<p->P125;++n)
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

void bedshear_probe::ini_location(lexer *p, ghostcell *pgc)
{
    int check;

    for(n=0;n<p->P125;++n)
    {
    iloc[n] = p->posc_i(p->P125_x[n]); 
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n] = p->posc_j(p->P125_y[n]); 

    check=ij_boundcheck(p,iloc[n],jloc[n],0);

	
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
