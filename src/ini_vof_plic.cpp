

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
    LOOP
    {
        if(p->pos_z()+0.5*p->DZN[KP]<p->F60)
            a->vof(i,j,k)=1.0;
        else if(p->pos_z()-0.5*p->DZN[KP]>p->F60)
        {
            a->vof(i,j,k)=0.0;
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
if(p->F100>0)
{
    double FSpos_pp,FSpos_pm,FSpos_mp,FSpos_mm,FSpos_c;
    double C_p,C_m;
    double No_z,No_y,No_x,St_loc,R_0,nsum,R_m;
    St_loc=p->F101_s;
    No_x=p->F101_nx;
    No_y=p->F101_ny;
    No_z=p->F101_nz;
    nsum=sqrt(No_x*No_x+No_y*No_y+No_z*No_z);
    No_x=No_x/nsum;
    No_y=No_y/nsum;
    No_z=No_z/nsum;
    if(p->F100==1)
    {
        R_0=No_x*St_loc;
        LOOP
        {
            FSpos_c=St_loc-No_y*p->pos_y()-No_z*p->pos_z();
            FSpos_pp=St_loc-No_y*(p->pos_y()+0.5*p->DYN[JP])-No_z*(p->pos_z()+0.5*p->DZN[KP]);
            FSpos_pm=St_loc-No_y*(p->pos_y()+0.5*p->DYN[JP])-No_z*(p->pos_z()-0.5*p->DZN[KP]);
            FSpos_mp=St_loc-No_y*(p->pos_y()-0.5*p->DYN[JP])-No_z*(p->pos_z()+0.5*p->DZN[KP]);
            FSpos_mm=St_loc-No_y*(p->pos_y()-0.5*p->DYN[JP])-No_z*(p->pos_z()-0.5*p->DZN[KP]);
            C_p=p->pos_x()+0.5*p->DXN[IP];
            C_m=p->pos_x()-0.5*p->DXN[IP];
        
            if((C_m<=FSpos_c&&FSpos_c<=C_p)||(C_m<=FSpos_pp&&FSpos_pp<=C_p)||(C_m<=FSpos_pm&&FSpos_pm<=C_p)||(C_m<=FSpos_mp&&FSpos_mp<=C_p)||(C_m<=FSpos_mm&&FSpos_mm<=C_m))
            {
                R_m=R_0-No_x*p->pos_x()-No_y*p->pos_y()-No_z*p->pos_z();
                a->vof(i,j,k)=VforPLIC(No_x,No_y,No_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],R_m);
            }
            else
            {
                if(p->pos_x()<=FSpos_c)
                {
                    if(No_x>=0.0)
                        a->vof(i,j,k)=1.0;
                    else
                        a->vof(i,j,k)=0.0;
                }
                else
                {
                    if(No_x>=0.0)
                        a->vof(i,j,k)=0.0;
                    else
                        a->vof(i,j,k)=1.0;
                }
            }
        }
    }
    else if(p->F100==2)
    {
        R_0=No_y*St_loc;
        LOOP
        {
            FSpos_c=St_loc-No_x*p->pos_x()-No_z*p->pos_z();
            FSpos_pp=St_loc-No_x*(p->pos_x()+0.5*p->DXN[IP])-No_z*(p->pos_z()+0.5*p->DZN[KP]);
            FSpos_mp=St_loc-No_x*(p->pos_x()-0.5*p->DXN[IP])-No_z*(p->pos_z()+0.5*p->DZN[KP]);
            FSpos_pm=St_loc-No_x*(p->pos_x()+0.5*p->DXN[IP])-No_z*(p->pos_z()-0.5*p->DZN[KP]);
            FSpos_mm=St_loc-No_x*(p->pos_x()-0.5*p->DXN[IP])-No_z*(p->pos_z()-0.5*p->DZN[KP]);
            C_p=p->pos_y()+0.5*p->DYN[JP];
            C_m=p->pos_y()-0.5*p->DYN[JP];
            
            if((C_m<=FSpos_c&&FSpos_c<=C_p)||(C_m<=FSpos_pp&&FSpos_pp<=C_p)||(C_m<=FSpos_pm&&FSpos_pm<=C_p)||(C_m<=FSpos_mp&&FSpos_mp<=C_p)||(C_m<=FSpos_mm&&FSpos_mm<=C_m))
            {
                R_m=R_0-No_x*p->pos_x()-No_y*p->pos_y()-No_z*p->pos_z();
                a->vof(i,j,k)=VforPLIC(No_x,No_y,No_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],R_m);
            }
            else
            {
                if(p->pos_y()<=FSpos_c)
                {
                    if(No_y>=0.0)
                        a->vof(i,j,k)=1.0;
                    else
                        a->vof(i,j,k)=0.0;
                }
                else
                {
                    if(No_y>=0.0)
                        a->vof(i,j,k)=0.0;
                    else
                        a->vof(i,j,k)=1.0;
                }
            }
        }
    }
    else
    {
        R_0=No_z*St_loc;
        LOOP
        {
            FSpos_c=St_loc-No_x*p->pos_x()-No_y*p->pos_y();
            FSpos_pp=St_loc-No_x*(p->pos_x()+0.5*p->DXN[IP])-No_y*(p->pos_y()+0.5*p->DYN[JP]);
            FSpos_mp=St_loc-No_x*(p->pos_x()-0.5*p->DXN[IP])-No_y*(p->pos_y()+0.5*p->DYN[JP]);
            FSpos_pm=St_loc-No_x*(p->pos_x()+0.5*p->DXN[IP])-No_y*(p->pos_y()-0.5*p->DYN[JP]);
            FSpos_mm=St_loc-No_x*(p->pos_x()-0.5*p->DXN[IP])-No_y*(p->pos_y()-0.5*p->DYN[JP]);
            C_p=p->pos_z()+0.5*p->DZN[KP];
            C_m=p->pos_z()-0.5*p->DZN[KP];
            if((C_m<=FSpos_c&&FSpos_c<=C_p)||(C_m<=FSpos_pp&&FSpos_pp<=C_p)||(C_m<=FSpos_pm&&FSpos_pm<=C_p)||(C_m<=FSpos_mp&&FSpos_mp<=C_p)||(C_m<=FSpos_mm&&FSpos_mm<=C_m))
            {
                R_m=R_0-No_x*p->pos_x()-No_y*p->pos_y()-No_z*p->pos_z();
                a->vof(i,j,k)=VforPLIC(No_x,No_y,No_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],R_m);
            }
            else
            {
                if(p->pos_z()<=FSpos_c)
                {
                    if(No_z>=0.0)
                        a->vof(i,j,k)=1.0;
                    else
                        a->vof(i,j,k)=0.0;
                }
                else
                {
                    if(No_z>=0.0)
                        a->vof(i,j,k)=0.0;
                    else
                        a->vof(i,j,k)=1.0;
                }
            }
        }
    }
    
}
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

double initialize::VforPLIC(double n_a, double n_b, double n_c, double d_a, double d_b, double d_c, double r0)
{
    double n_1, n_2, n_3, d_1, d_2, d_3, V, V0, vecsum, r;
    n_a=fabs(n_a);
    n_b=fabs(n_b);
    n_c=fabs(n_c);
    vecsum=sqrt(n_a*n_a+n_b*n_b+n_c*n_c);
    n_a=n_a/vecsum;
    n_b=n_b/vecsum;
    n_c=n_c/vecsum;
    
    if(n_b*d_b>=n_a*d_a)
    {
        n_1=n_a;
        d_1=d_a;
        n_2=n_b;
        d_2=d_b;
    }
    else
    {
        n_1=n_b;
        d_1=d_b;
        n_2=n_a;
        d_2=d_a;
    }
    
    if(n_c*d_c>=n_2*d_2)
    {
        n_3=n_c;
        d_3=d_c;
    }
    else if(n_c*d_c < n_1*d_1)
    {
        n_3=n_2;
        d_3=d_2;
        n_2=n_1;
        d_2=d_1;
        n_1=n_c;
        d_1=d_c;
    }
    else
    {
        n_3=n_2;
        d_3=d_2;
        n_2=n_c;
        d_2=d_c;
    }
    
    
    r=0.5*(n_1*d_1+n_2*d_2+n_3*d_3)-fabs(r0);
    if(r<=0.0) //case 0
    {
        V=0.0;
    }
    else if((min(n_1*d_1+n_2*d_2,n_3*d_3)<=r) && (r<=n_3*d_3)) //case 5
    {
        V=(r-0.5*(n_1*d_1+n_2*d_2))/(n_3*d_3);
    }
    else if(r<n_1*d_1) //case 1
    {
        V=(r*r*r)/(6.0*n_1*d_1*n_2*d_2*n_3*d_3);
    }
    else if(r<=n_2*d_2) //case 2
    {
        V=(3.0*r*(r-n_1*d_1)+n_1*n_1*d_1*d_1)/(6.0*n_2*d_2*n_3*d_3);
    }
    else //case 3&4
    {
        V=  (   r*r*r
                -(r-n_1*d_1)*(r-n_1*d_1)*(r-n_1*d_1)
                -(r-n_2*d_2)*(r-n_2*d_2)*(r-n_2*d_2)
                -fdim(r,n_3*d_3)*fdim(r,n_3*d_3)*fdim(r,n_3*d_3)    )
            /(6*n_1*d_1*n_2*d_2*n_3*d_3);
    }
    if(V!=V)
        cout<<"V in VOlcal NAN"<<endl;
    if(r0>=0)
    {
        V0=(0.5-V)+0.5;
    }
    else
    {
        V0=-(0.5-V)+0.5;
    }
    
    if(V0<0.0)
        cout<<"neg VO output"<<endl;
    if(V0>1.0)
        cout<<"too hight V0 output"<<endl;
        
    return V0;
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
