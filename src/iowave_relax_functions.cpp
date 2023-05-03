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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double iowave::rb1_ext(lexer *p, int var)
{
    double x0,y0,denom,r;
	double dist=1.0e20;
    int test1,test2,test_all;    
    int count;
    
    if(var==1)
    {
    x0 = p->pos1_x();
    y0 = p->pos1_y();
    }
    
    if(var==2)
    {
    x0 = p->pos2_x();
    y0 = p->pos2_y();
    }
    
    if(var==3||var==4)
    {
    x0 = p->pos_x();
    y0 = p->pos_y();
    }
    
    test_all=0;
    count=0;
    r = 0.0;
    
    for(int qn=0;qn<p->B108;++qn)
    {
    test1=0;
    test2=0;
    
    test1=intriangle(p,G1[qn][0],G1[qn][1],G3[qn][0],G3[qn][1],G2[qn][0],G2[qn][1],x0,y0);
    test2=intriangle(p,G3[qn][0],G3[qn][1],G4[qn][0],G4[qn][1],G2[qn][0],G2[qn][1],x0,y0);

        if(test1==1||test2==1)
        {
        test_all=1;
        
        // x dist
        denom = sqrt(pow(Ge[qn][1]-Gs[qn][1],2.0) + pow(Ge[qn][0]-Gs[qn][0],2.0));
        denom = denom>1.0e-20?denom:1.0e20;
        
        x = MIN(fabs((Ge[qn][1]-Gs[qn][1])*x0 - (Ge[qn][0]-Gs[qn][0])*y0 
                  + Ge[qn][0]*Gs[qn][1] - Ge[qn][1]*Gs[qn][0])/denom,dist);
        
        // relax
        dist1 = p->B108_d[qn]; 
        
        x=1.0-x/dist1;
        x=MAX(x,0.0);
        
        r += 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

        ++count;
        }
    }  
    
    if(test_all==0)
    r=1.0;
     
    if(test_all==1)
    r/=double(count);

    
      
    return r;
}

double iowave::rb3_ext(lexer *p, int var)
{
    double x0,y0,denom,r;
	double dist=1.0e20;
    int test1,test2,test_all;    
    int count;
    
    if(var==1)
    {
    x0 = p->pos1_x();
    y0 = p->pos1_y();
    }
    
    if(var==2)
    {
    x0 = p->pos2_x();
    y0 = p->pos2_y();
    }
    
    if(var==3||var==4)
    {
    x0 = p->pos_x();
    y0 = p->pos_y();
    }
    
    test_all=0;
    count=0;
    r = 0.0;
    
    for(int qn=0;qn<p->B107;++qn)
    {
    test1=0;
    test2=0;
    
    test1=intriangle(p,B1[qn][0],B1[qn][1],B3[qn][0],B3[qn][1],B2[qn][0],B2[qn][1],x0,y0);
    test2=intriangle(p,B3[qn][0],B3[qn][1],B4[qn][0],B4[qn][1],B2[qn][0],B2[qn][1],x0,y0);

        if(test1==1||test2==1)
        {
        test_all=1;
        
        // x dist
        denom = sqrt(pow(Be[qn][1]-Bs[qn][1],2.0) + pow(Be[qn][0]-Bs[qn][0],2.0));
        denom = denom>1.0e-20?denom:1.0e20;
        
        x = MIN(fabs((Be[qn][1]-Bs[qn][1])*x0 - (Be[qn][0]-Bs[qn][0])*y0 
                  + Be[qn][0]*Bs[qn][1] - Be[qn][1]*Bs[qn][0])/denom,dist);
        
        // relax
        
        dist2 = p->B107_d[qn]; 
        
        x=(dist2-fabs(x))/(dist2*dist2_fac);
        x=MAX(x,0.0);
        
        r += 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);
        ++count;
        }
        
    }
    
    
    if(test_all==0)
    r=1.0;
    
    if(test_all==1)
    r/=double(count);
	
	return r;
}

double iowave::rb1(lexer *p, double x)
{
    double r=0.0;

    x=1.0-x/dist1;
    x=MAX(x,0.0);
    
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);
      
    return r;
}

double iowave::rb3(lexer *p, double x)
{
    double r=0.0;

    x=(dist2-fabs(x))/(dist2*dist2_fac);
    x=MAX(x,0.0);
    
    
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);
	
	return r;
}

double iowave::ramp(lexer *p)
{
    double f=1.0;

    if(p->B101==1 && p->simtime<p->B102*p->wT)
    {
    f = p->simtime/(p->B102*p->wT) - (1.0/PI)*sin(PI*(p->simtime/(p->B102*p->wT)));
    }
    
    if(p->B101==2 && p->simtime<p->B102)
    {
    f = p->simtime/(p->B102) - (1.0/PI)*sin(PI*(p->simtime/(p->B102)));
    }

    return f;
}
