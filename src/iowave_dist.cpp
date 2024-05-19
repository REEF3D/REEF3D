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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double iowave::xgen(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos_x();
	y1 = p->pos_y();
	
	dist = fabs(y1 - tan_alpha*x1 + tan_alpha*x0 - y0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::xgen1(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos1_x();
	y1 = p->pos1_y();
	
	dist = fabs(y1 - tan_alpha*x1 + tan_alpha*x0 - y0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::xgen2(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos2_x();
	y1 = p->pos2_y();
	
	dist = fabs(y1 - tan_alpha*x1 + tan_alpha*x0 - y0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::ygen(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos_x();
	y1 = p->pos_y();
	
	dist = fabs(x1 - tan_alpha*y1 + tan_alpha*y0 - x0)/sqrt(pow(tan_alpha,2.0)+1.0);
    
	return dist;
}

double iowave::ygen1(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos1_x();
	y1 = p->pos1_y();
	
	dist = fabs(x1 - tan_alpha*y1 + tan_alpha*y0 - x0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::ygen2(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos2_x();
	y1 = p->pos2_y();
	
	dist = fabs(x1 - tan_alpha*y1 + tan_alpha*y0 - x0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::distgen(lexer *p)
{
    double x0,y0,denom;
	double dist=1.0e20;
    int test1,test2;    
    
    x0 = p->pos_x();
    y0 = p->pos_y();
    
    for(int qn=0;qn<p->B108;++qn)
    {
    test1=0;
    test2=0;
    
    test1=intriangle(p,G1[qn][0],G1[qn][1],G3[qn][0],G3[qn][1],G2[qn][0],G2[qn][1],x0,y0);
    test2=intriangle(p,G3[qn][0],G3[qn][1],G4[qn][0],G4[qn][1],G2[qn][0],G2[qn][1],x0,y0);

        if(test1==1||test2==1)
        {
        denom = sqrt(pow(Ge[qn][1]-Gs[qn][1],2.0) + pow(Ge[qn][0]-Gs[qn][0],2.0));
        denom = denom>1.0e-20?denom:1.0e20;
        
        dist = MIN(fabs((Ge[qn][1]-Gs[qn][1])*x0 - (Ge[qn][0]-Gs[qn][0])*y0 
                  + Ge[qn][0]*Gs[qn][1] - Ge[qn][1]*Gs[qn][0])/denom,dist);
        
        }
    }
    
	return dist;
}

double iowave::distbeach(lexer *p)
{
    double x0,y0,denom;
	double dist=1.0e20;
    int test1,test2;    
    
    x0 = p->pos_x();
    y0 = p->pos_y();
    
    for(int qn=0;qn<p->B107;++qn)
    {
    test1=0;
    test2=0;
    
    test1=intriangle(p,B1[qn][0],B1[qn][1],B3[qn][0],B3[qn][1],B2[qn][0],B2[qn][1],x0,y0);
    test2=intriangle(p,B3[qn][0],B3[qn][1],B4[qn][0],B4[qn][1],B2[qn][0],B2[qn][1],x0,y0);

        if(test1==1||test2==1)
        {
        denom = sqrt(pow(Be[qn][1]-Bs[qn][1],2.0) + pow(Be[qn][0]-Bs[qn][0],2.0));
        denom = denom>1.0e-20?denom:1.0e20;
        
        dist = MIN(fabs((Be[qn][1]-Bs[qn][1])*x0 - (Be[qn][0]-Bs[qn][0])*y0 
                  + Be[qn][0]*Bs[qn][1] - Be[qn][1]*Bs[qn][0])/denom,dist);
        
        }
    }
    
	return dist;
}

