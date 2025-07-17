

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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::inivofPLIC(fdm*a, lexer* p, ghostcell* pgc)
{

double dx=p->DXM;
double r;
double vofdiff, xdiff;

    LOOP
    {
	a->vof(i,j,k)=0.0;
    a->nZ(i,j,k)=1E06;
    a->nY(i,j,k)=1E06;
    a->nZ(i,j,k)=1E06;
    a->Alpha(i,j,k)=1E06;
    
    }
    
    pgc->start4(p,a->vof,1);

    if(p->F54<1E06 || p->F55 <1E06 || p->F56<1E06)
    {
        LOOP
        {
            if(p->XN[IP]>=p->F51 && p->XN[IP]<p->F54
            && p->YN[JP]>=p->F52 && p->YN[JP]<p->F55
            && p->ZN[KP]>=p->F53 && p->ZN[KP]<p->F56)
                a->vof(i,j,k)=1.0;
        }
    }

if(p->F57_1>1E-20||p->F57_2>1E-20||p->F57_3>1E-20||p->F57_4>1E-20)
{
	LOOP
	if(p->F57_1*p->XP[IP]+ p->F57_2*p->YP[JP]+ p->F57_3*p->ZP[KP] < p->F57_4)
        a->vof(i,j,k)=1.0;
}

if(p->F58_4>1E-20)
{
    LOOP
        {
        r = sqrt( pow(p->XP[IP]-p->F58_1,2.0)+pow(p->YP[JP]-p->F58_2,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));

        if(r<=p->F58_4)
        a->vof(i,j,k)=1.0;
        }
}

if(p->F60>-1.0e20)
{   p->phimean=p->F60;
    cout<<"a"<<endl;
    LOOP
    {
        if(p->pos_z()+0.5*p->DZN[KP]<p->F60)
            a->vof(i,j,k)=1.0;
        else if(p->pos_z()-0.5*p->DZN[KP]>p->F60)
        {
            a->vof(i,j,k)=0.0;
            cout<<"b"<<endl;
        }
        else
        {
            a->vof(i,j,k)=(p->F60-(p->pos_z()-0.5*p->DZN[KP]))/p->DZN[KP];
        }
    }
}

  /*  if((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20&& p->F63>-1.0e-20  )
    {
        vofdiff=p->F62-p->phimean;
        xdiff=p->xcoormax-p->F63;

        LOOP
        if(p->pos_x() > p->F63)
        a->vof(i,j,k)=(vofdiff/xdiff)*(p->pos_x()-p->F63) + p->phimean    - p->pos_z() ;
    }*/

	double H=0.0;

	LOOP
	{
		H=a->vof(i,j,k);

		H=MAX(H,0.0);
		H=MIN(H,1.0);

		a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}
    
    //LOOP
    //a->phi(i,j,k) = a->vof(i,j,k);

	pgc->start4(p,a->vof,50);
   // pgc->start4(p,a->phi,50);
	pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);
}



/*

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::inivofPLIC(fdm*a, lexer* p, ghostcell* pgc)
{

    double dx=p->DXM;
    double r;
    double vofdiff, xdiff;
    
    p->phimean=p->F56;

    LOOP
	a->vof(i,j,k)=0.0;
    
	LOOP
	if 
    (
        double(i)*dx + p->originx >= p->F51 && double(i)*dx + p->originx < p->F54
	 && double(j)*dx + p->originy >= p->F52 && double(j)*dx + p->originy < p->F55
    )
    {
        double value;    
        LOOP
        {
            if (p->phimean >= p->pos_z() + p->DZN[KP]/2.0)
            {
                value = 1.0;
            }
            else if (p->phimean <= p->pos_z() - p->DZN[KP]/2.0)
            {
                value = 0.0;
            }
            else 
            {
                value = (p->phimean - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
            }
            
            a->vof(i,j,k) = value;
        }
    }

    if (p->F57_1 > 0 || p->F57_2 > 0 || p->F57_3 > 0 || p->F57_4 > 0)
    {
        LOOP
        if 
        (
            p->F57_1*((double(i) + 0.5)*dx + p->originx) 
          + p->F57_2*((double(j) + 0.5)*dx + p->originy)
          + p->F57_3*((double(k) + 0.5)*dx + p->originz) 
          < p->F57_4
        )
        {
            a->vof(i,j,k)=1.0;
        }
    }

    if(p->F58_4>0.0)
    {
        p->F58_1 -= p->originx;
        p->F58_2 -= p->originy;
        p->F58_3 -= p->originz;

        LOOP
        {
            r = 
                sqrt
                ( 
                    pow((double(i) + 0.5)*dx - p->F58_1, 2.0)
                  + pow((double(j) + 0.5)*dx - p->F58_2, 2.0)
                  + pow((double(k) + 0.5)*dx - p->F58_3, 2.0)
                );
            
            if(r<=p->F58_4)
            a->vof(i,j,k)=1.0;
        }
    }

    if (p->F60 > -1.0e20)
    {        
        p->phimean=p->F60;
        
        double value;
        
        LOOP
        {
            if (p->phimean >= p->pos_z() + p->DZN[KP]/2.0)
            {
                value = 1.0;
            }
            else if (p->phimean <= p->pos_z() - p->DZN[KP]/2.0)
            {
                value = 0.0;
            }
            else 
            {
                value = (p->phimean - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
            }
            
            a->vof(i,j,k) = value;
        }
    }


    if ((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20 && p->F63>-1.0e-20)
    {
        vofdiff = p->F62-p->phimean;
        xdiff = p->xcoormax-p->F63;

        LOOP
        if (p->pos_x() > p->F63)
        a->vof(i,j,k) = (vofdiff/xdiff)*(p->pos_x()-p->F63) + p->phimean - p->pos_z();
    }

	double H=0.0;

	LOOP
	{
		H = a->vof(i,j,k);

		H = MAX(H, 0.0);
		H = MIN(H, 1.0);

		a->ro(i,j,k) = p->W1*H + p->W3*(1.0 - H);
		a->visc(i,j,k) = p->W2*H + p->W4*(1.0 - H);
	}
	pgc->start4(p,a->vof,50);
	pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);
*/
/*
    //- Initialise distance function at start of simulation
    p->F40 = 23;
    if(p->F70 > 0 || p->F71 > 0 || p->F72 > 0)
    {
        iniphi_box(p, a, pgc);
    }
    else
    {
        iniphi(a, p, pgc);
    } 

	LOOP
	{
		a->test(i,j,k) = a->vof(i,j,k);
	} 
	pgc->start4(p,a->test,50);
}
*/
