

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Architect: Hans Bihs
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
    a->nX(i,j,k)=1E06;
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

if(p->F100>0)
{
    // plane with the normal F 101 nx ny nz (pointing out of the water) through the point at F 101 e on
    // the x (F 100 1), y (F 100 2) or z axis (else). VforPLIC gives the exact fraction of every cell,
    // also 1 or 0 for the cells the plane does not cut
    double No_z,No_y,No_x,St_loc,R_0,nsum,R_m;
    St_loc=p->F101_e;
    No_x=p->F101_nx;
    No_y=p->F101_ny;
    No_z=p->F101_nz;
    nsum=sqrt(No_x*No_x+No_y*No_y+No_z*No_z);
    No_x=No_x/nsum;
    No_y=No_y/nsum;
    No_z=No_z/nsum;
    
    if(p->F100==1)
        R_0=No_x*St_loc;
    else if(p->F100==2)
        R_0=No_y*St_loc;
    else
        R_0=No_z*St_loc;
        
    LOOP
    {
        R_m=R_0-No_x*p->pos_x()-No_y*p->pos_y()-No_z*p->pos_z();
        a->vof(i,j,k)=VforPLIC(No_x,No_y,No_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],R_m);
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

	pgc->start4(p,a->vof,50);
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
