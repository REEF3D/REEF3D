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

#include"sflow_print_bed.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_print_bed::sflow_print_bed(lexer *p, fdm2D* b)
{

	gauge_num = p->P121;
	x = p->P121_x;
	y = p->P121_y;

	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SFLOW_BED",0777);
	
    if(p->mpirank==0 && p->P121>0)
    {
    // open file
	if(p->P14==0)
    bedout.open("REEF3D-SFLOW-BED-HG.dat");
	
	if(p->P14==1)
	bedout.open("./REEF3D_SFLOW_BED/REEF3D-SFLOW-BED-HG.dat");

    bedout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    bedout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    bedout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    bedout<<endl<<endl;

    bedout<<"time";
    for(n=0;n<gauge_num;++n)
    bedout<<"\t P"<<n+1;

    bedout<<endl<<endl;
    }
	
	//-------------------
	
	
	p->Iarray(iloc,gauge_num);
	p->Iarray(jloc,gauge_num);
	p->Iarray(flag,gauge_num);
	p->Darray(bed,gauge_num);

    ini_location(p,b);
}

sflow_print_bed::~sflow_print_bed()
{
    bedout.close();
}

void sflow_print_bed::height_gauge(lexer *p, fdm2D *b, ghostcell *pgc, slice &f)
{
    double zval=0.0;

    for(n=0;n<gauge_num;++n)
    bed[n]=-1.0e20;

	
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
			bed[n] = f(i,j);

    }
	
    for(n=0;n<gauge_num;++n)
    bed[n]=pgc->globalmax(bed[n]);

    // write to file
    if(p->mpirank==0)
    {
    bedout<<setprecision(9)<<p->simtime<<"\t";
    for(n=0;n<gauge_num;++n)
    bedout<<setprecision(9)<<bed[n]<<"  \t  ";
    bedout<<endl;
    }
}

void sflow_print_bed::ini_location(lexer *p, fdm2D *b)
{
    for(n=0;n<gauge_num;++n)
    {
    iloc[n]=conv((x[n]-p->originx)/p->DXM);
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n]=conv((y[n]-p->originy)/p->DXM);

    if(iloc[n]>=0 && iloc[n]<p->knox)
    if(jloc[n]>=0 && jloc[n]<p->knoy)
    flag[n]=1;
    }
}

int sflow_print_bed::conv(double a)
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
