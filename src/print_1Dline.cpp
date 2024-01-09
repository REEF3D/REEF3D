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

#include"print_1Dline.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

print_1Dline::print_1Dline(lexer *p, fdm* a, ghostcell *pgc)
{	
	p->Iarray(iloc,p->P51);
	p->Iarray(jloc,p->P51);
	p->Iarray(flag,p->P51);
	p->Darray(wsf,p->P51);

    for(n=0;n<p->P51;++n)
    {
    iloc[n]=0;
    jloc[n]=0;
    flag[n]=0;
    wsf[n]=0.0;
    }

    if(p->mpirank==0)
    {
    // open file
    wsfout.open("REEF3D_WSF-HG.dat");

    wsfout<<"number of gauges:  "<<p->P51<<endl<<endl;
    wsfout<<"x_coord     ycoord"<<endl;
    for(n=0;n<p->P51;++n)
    wsfout<<n+1<<"\t "<<p->P51_x[n]<<"\t "<<p->P51_y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<p->P51;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }

    ini_location(p,a,pgc);
}

print_1Dline::~print_1Dline()
{
    wsfout.close();
}

void print_1Dline::height_gauge(lexer *p, fdm *a, ghostcell *pgc)
{
    double zval=0.0;

    for(n=0;n<p->P51;++n)
    wsf[n]=-1.0e20;


    for(n=0;n<p->P51;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];

        KLOOP
        PCHECK
        {
            if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
            wsf[n]=MAX(wsf[n],-(a->phi(i,j,k)*p->DXM)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->ZP[KP]);
        }
    }

    for(n=0;n<p->P51;++n)
    wsf[n]=pgc->globalmax(wsf[n]);

    // write to file
    if(p->mpirank==0)
    {
    wsfout<<p->simtime<<"\t";
    for(n=0;n<p->P51;++n)
    wsfout<<wsf[n]<<"  \t  ";
    wsfout<<endl;
    }
}

void print_1Dline::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void print_1Dline::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check;

    for(n=0;n<p->P51;++n)
    {
    iloc[n]=conv((p->P51_x[n]-p->originx)/p->DXM);
    jloc[n]=conv((p->P51_y[n]-p->originy)/p->DXM);

    check=ij_boundcheck(p,a,iloc[n],jloc[n],0);

    if(check==1)
    flag[n]=1;
    }
}


int print_1Dline::conv(double a)
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


