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

#include"flowfile_in.h"
#include"lexer.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

double flowfile_in::ccipol4(lexer *p, double**F, double xp, double yp, double zp)
{
    int ii,jj,kk;
    double wa,wb,wc,value;
    
    ii=i;
    jj=j;
    kk=k;
    
    xp -= p->originx;
    yp -= p->originy;
    zp -= p->originz;

    i=0;
    j=0;
    
    k=int((zp-zs)/deltax-0.5);
	

    k=MAX(k,0);
    k=MIN(k,Nk);
    


    wa=((double(i) + 1.5)-xp/deltax);
    wb=((double(j) + 1.5)-yp/deltax);
    wc=((double(k) + 1.5)-zp/deltax);
	

    value =  lint4(F,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double flowfile_in::lint4(double **F, int& i,int& j, int& k, double wa, double wb, double wc)
{
    double y1,y2,value;

    /*x1 = wa*f(i,j,k)   + (1.0-wa)*f(i+1,j,k);
    x2 = wa*f(i,j+1,k) + (1.0-wa)*f(i+1,j+1,k);

    x3 = wa*f(i,j,k+1)   + (1.0-wa)*f(i+1,j,k+1);
    x4 = wa*f(i,j+1,k+1) + (1.0-wa)*f(i+1,j+1,k+1);

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;*/
    
    y1 = F[0][k];
    y2 = F[0][k+1];

    value = wc*y1 +(1.0-wc)*y2;
	
    
    return value;
}
