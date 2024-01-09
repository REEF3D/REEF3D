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

#include"interpolation.h"
#include"field.h"
#include"lexer.h"
#include"fdm.h"

double interpolation::lintsl1(slice &f, int& i,int& j, double wa, double wb)
{
    v1=v2=v3=v4=0.0;
    c1=c2=c3=c4=0;

pip=4;
    if(p->flagslice1[IJ]>0)
    {
    v1=f(i,j);
    c1=1;
    }
    
    if(p->flagslice1[IJp1]>0)
    {
    v2=f(i,j+1);
    c2=1;
    }
    
    if(p->flagslice1[Ip2J]>0)
    {
    v3=f(i+1,j);
    c3=1;
    }
    
    if(p->flagslice1[Ip1Jp1]>0)
    {
    v4=f(i+1,j+1);
    c4=1;
    }
pip=0;
    
    // x1
    if(c1==1 && c3==1)
    x1 = wa*v1 + (1.0-wa)*v3;
    
    if(c1==1 && c3==0)
    x1 = v1;
    
    if(c1==0 && c3==1)
    x1 = v3;
    
    
    // x2
    if(c2==1 && c4==1)
    x2 = wa*v2 + (1.0-wa)*v4;
    
    if(c2==1 && c4==0)
    x2 = v2;
    
    if(c2==0 && c4==1)
    x2 = v4;
    
    if((c1==0 && c3==0) && (c2==1 || c4==1))
    wb=0.0;
    
    if((c2==0 && c4==0) && (c1==1 || c3==1))
    wb=1.0;
    
    if(c2==0 && c4==0 && c1==1 && c3==1)
    {
    x1=x2=0.0;
    }

    value = wb*x1 + (1.0-wb)*x2;
    
 return value;
}

double interpolation::lintsl2(slice &f, int& i,int& j, double wa, double wb)
{
        v1=v2=v3=v4=0.0;
    c1=c2=c3=c4=0;

pip=4;
    if(p->flagslice2[IJ]>0)
    {
    v1=f(i,j);
    c1=1;
    }
    
    if(p->flagslice2[IJp1]>0)
    {
    v2=f(i,j+1);
    c2=1;
    }
    
    if(p->flagslice2[Ip2J]>0)
    {
    v3=f(i+1,j);
    c3=1;
    }
    
    if(p->flagslice2[Ip1Jp1]>0)
    {
    v4=f(i+1,j+1);
    c4=1;
    }
pip=0;
    
    // x1
    if(c1==1 && c3==1)
    x1 = wa*v1 + (1.0-wa)*v3;
    
    if(c1==1 && c3==0)
    x1 = v1;
    
    if(c1==0 && c3==1)
    x1 = v3;
    
    
    // x2
    if(c2==1 && c4==1)
    x2 = wa*v2 + (1.0-wa)*v4;
    
    if(c2==1 && c4==0)
    x2 = v2;
    
    if(c2==0 && c4==1)
    x2 = v4;
    
    if((c1==0 && c3==0) && (c2==1 || c4==1))
    wb=0.0;
    
    if((c2==0 && c4==0) && (c1==1 || c3==1))
    wb=1.0;
    
    if(c2==0 && c4==0 && c1==1 && c3==1)
    {
    x1=x2=0.0;
    }

    value = wb*x1 + (1.0-wb)*x2;
    
 return value;
}


double interpolation::lintsl4(slice& f, int& i,int& j, double wa, double wb)
{
    v1=v2=v3=v4=0.0;
    c1=c2=c3=c4=0;

pip=4;
    if(p->flagslice4[IJ]>0)
    {
    v1=f(i,j);
    c1=1;
    }
    
    if(p->flagslice4[IJp1]>0)
    {
    v2=f(i,j+1);
    c2=1;
    }
    
    if(p->flagslice4[Ip2J]>0)
    {
    v3=f(i+1,j);
    c3=1;
    }
    
    if(p->flagslice4[Ip1Jp1]>0)
    {
    v4=f(i+1,j+1);
    c4=1;
    }
pip=0;
    
    // x1
    if(c1==1 && c3==1)
    x1 = wa*v1 + (1.0-wa)*v3;
    
    if(c1==1 && c3==0)
    x1 = v1;
    
    if(c1==0 && c3==1)
    x1 = v3;
    
    
    // x2
    if(c2==1 && c4==1)
    x2 = wa*v2 + (1.0-wa)*v4;
    
    if(c2==1 && c4==0)
    x2 = v2;
    
    if(c2==0 && c4==1)
    x2 = v4;
    
    if((c1==0 && c3==0) && (c2==1 || c4==1))
    wb=0.0;
    
    if((c2==0 && c4==0) && (c1==1 || c3==1))
    wb=1.0;
    
    /*
    if(c2==0 && c4==0 && c1==1 && c3==1)
    {
    //x1=x2=0.0;
    }*/

    value = wb*x1 + (1.0-wb)*x2;
    
 return value;

}
