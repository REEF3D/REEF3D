/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"bedprobe_point.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

bedprobe_point::bedprobe_point(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    p->Iarray(iloc,p->P121);
	p->Iarray(jloc,p->P121);
	p->Iarray(flag,p->P121);
	p->Darray(wsf,p->P121);
	
	// Create Folder
	if(p->mpirank==0 && p->A10==2)
	mkdir("./REEF3D_SFLOW_Sediment",0777);
    
    if(p->mpirank==0 && p->A10==5)
	mkdir("./REEF3D_NHFLOW_Sediment",0777);
    
    if(p->mpirank==0 && p->A10==6)
	mkdir("./REEF3D_CFD_Sediment",0777);
	
    if(p->mpirank==0 && p->P121>0)
    {
    // open file
    if(p->A10==2)
	wsfout.open("./REEF3D_SFLOW_Sediment/REEF3D-SFLOW-Sediment-Point.dat");
    
    if(p->A10==5)
	wsfout.open("./REEF3D_NHFLOW_Sediment/REEF3D-NHFLOW-Sediment-Point.dat");
    
    if(p->A10==6)
	wsfout.open("./REEF3D_CFD_Sediment/REEF3D-CFD-Sediment-Point.dat");

    wsfout<<"number of gauges:  "<<p->P121<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<p->P121;++n)
    wsfout<<n+1<<"\t "<<p->P121_x[n]<<"\t "<<p->P121_y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<p->P121;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }

    ini_location(p,pgc,s);
}

bedprobe_point::~bedprobe_point()
{
    wsfout.close();
}

void bedprobe_point::bed_gauge(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double zval=0.0;

    for(n=0;n<p->P121;++n)
    wsf[n]=-1.0e20;

	
    for(n=0;n<p->P121;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
    wsf[n] = MAX(wsf[n],s->bedzh(i,j));
    }
	
    for(n=0;n<p->P121;++n)
    wsf[n]=pgc->globalmax(wsf[n]);

    // write to file
    if(p->mpirank==0)
    {
    wsfout<<p->sedtime<<"\t ";
    for(n=0;n<p->P121;++n)
    wsfout<<wsf[n]<<"  \t  ";
    wsfout<<endl;
    }
}

void bedprobe_point::write(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
}

void bedprobe_point::ini_location(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    int check;

    for(n=0;n<p->P121;++n)
    {
    iloc[n] = p->posc_i(p->P121_x[n]);
    jloc[n] = p->posc_j(p->P121_y[n]);

    check=ij_boundcheck(p,iloc[n],jloc[n],0);

    if(check==1)
    flag[n]=1;
	
	cout<<p->mpirank<<" n: "<<n<<" x: "<<p->P121_x[n]<<" y: "<<p->P121_y[n]<<" iloc: "<<iloc[n]<<" jloc: "<<jloc[n]<<" n: "<<n<<" flag: "<<flag[n]<<endl;
    }
}


int bedprobe_point::conv(double a)
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
