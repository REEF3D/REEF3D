/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm2D.h"
#include"slice.h"
             
double ghostcell::gcsl_ipol1(lexer* p, fdm2D *b, slice &f)
{
    v1=v2=0;

    pip=4;
    if(p->flagslice1[IJ]>0)
    v1=f(i,j);
    if(p->flagslice1[IJp1]>0)
    v2=f(i,j+1);
    pip=0;

    value = 0.5*(v1+v2);

    return value;
}

double ghostcell::gcsl_ipol2(lexer* p, fdm2D *b, slice &f)
{
    v1=v2=0;

    pip=4;
    if(p->flagslice2[IJ]>0)
    v1=f(i,j);
    if(p->flagslice2[Ip1J]>0)
    v2=f(i+1,j);
    pip=0;

    value = 0.5*(v1+v2);

    return value;
}

double ghostcell::gcsl_ipol1a(lexer* p, fdm2D *b, slice &f)
{
    v1=v2=0;

    pip=4;
    v1=f(i,j);
    v2=f(i,j+1);
    pip=0;

    value = 0.5*(v1+v2);

    return value;
}

double ghostcell::gcsl_ipol2a(lexer* p, fdm2D *b, slice &f)
{
    v1=v2=0;

    pip=4;
    v1=f(i,j);
    v2=f(i+1,j);
    pip=0;

    value = 0.5*(v1+v2);

    return value;
}

double ghostcell::gcsl_ipol4(lexer* p, slice &f)
{
    v1=v2=v3=v4=0;

    pip=4;
    v1=f(i,j);

    v2=f(i+1,j);

    v3=f(i,j+1);

    v4=f(i+1,j+1);
    pip=0;

    value = 0.25*(v1+v2+v3+v4);
    
    return value;
}

double ghostcell::gcsl_ipol4eta(lexer* p, fdm2D *b, slice &f)
{
    double bedvalue;
    
    double wd_criterion=0.00005;
    
    if(p->A244==1)
    wd_criterion=p->A244_val;
    
    if(p->A245==1)
    wd_criterion=p->A245_val*p->dx;
    
    v1=v2=v3=v4=0;

    pip=4;
    v1=f(i,j);

    v2=f(i+1,j);

    v3=f(i,j+1);

    v4=f(i+1,j+1);
    pip=0;

    value = 0.25*(v1+v2+v3+v4);
    
    //bed
    pip=4;
    v1=b->bed(i,j);

    v2=b->bed(i+1,j);

    v3=b->bed(i,j+1);

    v4=b->bed(i+1,j+1);
    pip=0;

    bedvalue = 0.25*(v1+v2+v3+v4);
    
    if(value+p->wd>bedvalue)
    if(value+p->wd-bedvalue<wd_criterion)
    value=value-1.0*wd_criterion;
    
  

    return value;
}

double ghostcell::gcsl_ipolint(lexer* p, fdm2D *b, sliceint &f)
{
    v1=v2=v3=v4=0;

    pip=4;
    v1=double(f(i,j));

    v2=double(f(i+1,j));

    v3=double(f(i,j+1));

    v4=double(f(i+1,j+1));
    pip=0;

    value = 0.25*(v1+v2+v3+v4);

    return value;
}

double ghostcell::gcsl_ccipol4(lexer* p,slice& f, double xp, double yp)
{
    ii=i;
    jj=j;

    i=int((xp-p->originx)/dx-0.5);
		if((xp-p->originx)/dx-0.5<0.0)
		--i;
    j=int((yp-p->originy)/dx-0.5);
		if((yp-p->originy)/dx-0.5<0.0)
		--j;

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);


    value =  gcsl_lint4(p,f,i,j,wa,wb);

    i=ii;
    j=jj;

    return value;
}

double ghostcell::gcsl_lint4(lexer *p, slice& f, int& i,int& j, double wa, double wb)
{
    v1=v2=v3=v4=0.0;

pip=4;

    v1=f(i,j);

    v2=f(i,j+1);

    v3=f(i+1,j);

    v4=f(i+1,j+1);

    pip=0;
pip=0;

    x1 = wa*v1   + (1.0-wa)*v3;
    x2 = wa*v2   + (1.0-wa)*v4;

    value = wb*x1 +(1.0-wb)*x2;

pip=0;
 return value;

}