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

#include"bedprobe_point.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

bedprobe_point::bedprobe_point(lexer *p, fdm* a, ghostcell *pgc)
{
    p->Iarray(iloc,p->P121);
	p->Iarray(jloc,p->P121);
	p->Iarray(flag,p->P121);
	p->Darray(wsf,p->P121);
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_SedimentPoint",0777);
	
    if(p->mpirank==0 && p->P121>0)
    {
    // open file
	if(p->P14==0)
    wsfout.open("REEF3D-CFD-Sediment-Point.dat");
	
	if(p->P14==1)
	wsfout.open("./REEF3D_CFD_SedimentPoint/REEF3D-CFD-Sediment-Point.dat");

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

    ini_location(p,a,pgc);
}

bedprobe_point::~bedprobe_point()
{
    wsfout.close();
}

void bedprobe_point::bed_gauge(lexer *p, fdm *a, ghostcell *pgc)
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
	
	//cout<<p->mpirank<<" n: "<<n<<" flag: "<<flag[n]<<" iloc: "<<iloc[n]<<" jloc: "<<jloc[n]<<endl;

        KLOOP
        PBASECHECK
        {
            if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
            wsf[n]=MAX(wsf[n],-(a->topo(i,j,k)*p->DXM)/(a->topo(i,j,k+1)-a->topo(i,j,k)) + p->pos_z());
        }
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

void bedprobe_point::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void bedprobe_point::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check;

    for(n=0;n<p->P121;++n)
    {
    iloc[n]=conv((p->P121_x[n]-p->originx)/p->DXM);
    jloc[n]=conv((p->P121_y[n]-p->originy)/p->DXM);

    check=ij_boundcheck(p,a,iloc[n],jloc[n],0);

    if(check==1)
    flag[n]=1;
	
	//cout<<p->mpirank<<" n: "<<n<<" x: "<<p->P121_x[n]<<" y: "<<p->P121_y[n]<<" iloc: "<<iloc[n]<<" jloc: "<<jloc[n]<<" n: "<<n<<" flag: "<<flag[n]<<endl;
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
