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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::updateVOF(fdm* a, lexer* p, int sweep)
{
	if (sweep == 0)
	{
		LOOP
		{   
            if(a->vof(i,j,k)<=0.001)
                a->vof(i,j,k)=(vof2(i,j,k)+vof1w(i+1,j,k)+vof3w(i-1,j,k)+vof1s(i+1,j,k)+vof3s(i-1,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
            else if(a->vof(i,j,k)>=0.999)
            {
                if(vof1s(i+1,j,k)<=0.00001 && vof3s(i-1,j,k)<=0.000001)
                    a->vof(i,j,k)=vof2(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else if(vof1s(i+1,j,k)>0.00001 && vof3s(i-1,j,k)<=0.00001)
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q1(i+1,j,k))*p->DYN[JP]*p->DZN[KP]-vof1s(i+1,j,k)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else if(vof1s(i+1,j,k)<=0.000001 && vof3s(i-1,j,k)>0.00001)
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q2(i-1,j,k))*p->DYN[JP]*p->DZN[KP]-vof3s(i-1,j,k)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q1(i+1,j,k))*p->DYN[JP]*p->DZN[KP]-vof1s(i+1,j,k))-(abs(Q2(i-1,j,k))*p->DYN[JP]*p->DZN[KP]-vof3s(i-1,j,k)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
            else
                a->vof(i,j,k)=(vof2(i,j,k)+vof1s(i+1,j,k)+vof3s(i-1,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
        a->vof(i,j,k)=max(0.0,min(a->vof(i,j,k),1.0));
		}
	}
	else if (sweep == 1)
	{
		LOOP
		{   
            if(a->vof(i,j,k)<=0.001)
                a->vof(i,j,k)=(vof2(i,j,k)+vof1w(i,j+1,k)+vof3w(i,j-1,k)+vof1s(i,j+1,k)+vof3s(i,j-1,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
            else if(a->vof(i,j,k)>=0.999)
            {
                if(vof1s(i,j+1,k)<=1e-06 && vof3s(i,j-1,k)<=1e-06)
                    a->vof(i,j,k)=vof2(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else if(vof1s(i,j+1,k)>1e-06 && vof3s(i,j-1,k)<=1e-06)
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q1(i,j+1,k))*p->DXN[IP]*p->DZN[KP]-vof1s(i,j+1,k)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else if(vof1s(i,j+1,k)<=1e-06 && vof3s(i,j-1,k)>1e-06)
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q2(i,j-1,k))*p->DXN[IP]*p->DZN[KP]-vof3s(i,j-1,k)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q1(i,j+1,k))*p->DXN[IP]*p->DZN[KP]-vof1s(i,j+1,k))-(abs(Q2(i,j-1,k))*p->DXN[IP]*p->DZN[KP]-vof3s(i,j-1,k)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
            else
                a->vof(i,j,k)=(vof2(i,j,k)+vof1s(i,j+1,k)+vof3s(i,j+1,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
        a->vof(i,j,k)=max(0.0,min(a->vof(i,j,k),1.0));
		}
	}
	else
	{
        LOOP
        {
		if(a->vof(i,j,k)<=0.001)
                a->vof(i,j,k)=(vof2(i,j,k)+vof1w(i,j,k+1)+vof3w(i,j,k-1)+vof1s(i,j,k+1)+vof3s(i,j,k-1))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
            else if(a->vof(i,j,k)>=0.999)
            {
                if(vof1s(i,j,k+1)<=0.00001 && vof3s(i,j,k-1)<=0.00001)
                    a->vof(i,j,k)=vof2(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else if(vof1s(i,j,k+1)>0.00001 && vof3s(i,j,k-1)<=0.00001)
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q1(i,j,k+1))*p->DXN[IP]*p->DYN[JP]-vof1s(i,j,k+1)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else if(vof1s(i,j,k+1)<=0.00001 && vof3s(i,j,k-1)>0.00001)
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q2(i,j,k-1))*p->DXN[IP]*p->DYN[JP]-vof3s(i,j,k-1)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                    
                else
                    a->vof(i,j,k)=(vof2(i,j,k)-(abs(Q1(i,j,k+1))*p->DXN[IP]*p->DYN[JP]-vof1s(i,j,k+1))-(abs(Q2(i,j,k-1))*p->DXN[IP]*p->DYN[JP]-vof3s(i,j,k-1)))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
            else
                a->vof(i,j,k)=(vof2(i,j,k)+vof1s(i,j,k+1)+vof3s(i,j,k-1))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
        a->vof(i,j,k)=max(0.0,min(a->vof(i,j,k),1.0));
		}

	}
}


void VOF_PLIC::updateVolumeFraction
(
	fdm* a, 
	lexer* p,
	int sweep
)
{
    vof1w(i,j,k)=0.0;
    vof3w(i,j,k)=0.0;
	//- Calculate volume entering, leaving and staying in the cell
    if(sweep==0)
    {
    if(Q1(i,j,k)<0.0)
        vof1s(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,1);
    else
        vof1s(i,j,k)=0.0;
        
    if(Q2(i,j,k)>0.0)
        vof3s(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,3);
    else
        vof3s(i,j,k)=0.0;
    
    vof2(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,2);
    }
    
    if(sweep==1)
    {
    if(Q1(i,j,k)<-0.0)
        vof1s(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,1);
    else
        vof1s(i,j,k)=0.0;
        
    if(Q2(i,j,k)>0.0)
        vof3s(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,3);
    else
        vof3s(i,j,k)=0.0;
    
    vof2(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,2);
    }
    
    if(sweep==2)
    {
    if(Q1(i,j,k)<0.0)
        vof1s(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,1);
    else
        vof1s(i,j,k)=0.0;
        
    if(Q2(i,j,k)>0.0)
        vof3s(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,3);
    else
        vof3s(i,j,k)=0.0;
    
    vof2(i,j,k)=Volume_flow(p,nx(i,j,k),ny(i,j,k),nz(i,j,k),alpha(i,j,k),sweep,2);
    }
    
}

double VOF_PLIC::calcV
(
     double& M1,
     double& m2,
     double& m3,
     double& alpha,
	double r0,
	double dr0
)
{
	//- Move origin to r0 along axis as given in Gueyffier (25)
	// if Q1 < 0: vof1 starts at r0=Q1 and goes to  0 --> dr0 = -Q1
	// if Q2 > 0: vof3 starts at r0= 1 and goes to Q2 --> dr0 = Q2
	// vof2 starts at r0=0 or r0=Q1 (if Q1 > 0) and ends at 1, 1-Q1 (if Q1 > 0), 1-Q1-abs(Q2) (if Q2 < 0)
	
	double al = alpha - M1*r0;

	
	//- Reflect parallelepiped such that m_i > 0
	
	al += max(0.0, -M1*dr0) + max(0.0, -m2) + max(0.0, -m3);
	
	
	//- Normalise plane equation: n1*x + n2*y + n3*z = al
	
	double tmp = sqrt(fabs(M1)*dr0*fabs(M1)*dr0 + fabs(m2)*fabs(m2) + fabs(m3)*fabs(m3));
	if (tmp < 1e-10) return 0.0;	
	
	double n1 = fabs(M1)/(tmp + 1e-50);
	double n2 = fabs(m2)/(tmp + 1e-50);
	double n3 = fabs(m3)/(tmp + 1e-50);
	
	al = max(0.0, min(1.0, al/tmp));
	
	
	//- Limit alpha such that 0.0 < alpha < 0.5
	
	double al0 = min(al, 1.0 - al);
	
	
	//- Order b_i such that b1 < b2 < b3
	
	double b1 = min(n1*dr0, n2);
	double b3 = max(n1*dr0, n2);
	double b2 = n3;
	if (b2 < b1)
	{
		tmp = b1;
		b1 = b2;
		b2 = tmp;
	}
	else if (b2 > b3)
	{
		tmp = b3;
		b3 = b2;
		b2 = tmp;
	}
	double b12 = b1 + b2;
	double bm = min(b12, b3);
	
	
	//- Calculate new volume, assuming x0 = 0.0 and delta_x = 1.0 
	// Scardovelli p.233: alpha = al0, M1 = b1, m2 = m2, m3 = b3, m = bm, M12 = b12, m3 = b3, V = tmp
	
	double pr = max(6.0*b1*b2*b3, 1.0e-20);
	
	// if alpha < b1
	if (al0 < b1)
	{
		tmp = pow(al0, 3.0)/pr;
	}
	else if (al0 < b2)	// if b1 < alpha < b2
	{
		tmp = 0.5*al0*(al0 - b1)/(b2*b3 + 1e-20) + pow(b1, 3.0)/pr;
	}
	else if (al0 < bm)	// if b2 < alpha < bm
	{
		tmp = 
			(al0*al0*(3.0*b12 - al0) + b1*b1*(b1 - 3.0*al0) 
			+ b2*b2*(b2 - 3.0*al0))
			/pr;
	}
	else if (b12 < b3)	// if b12 < b3 --> bm = b12
	{
		tmp = (al0 - 0.5*bm)/(b3 + 1e-20);
	}
	else	// if b3 < b12 --> bm = b3
	{
		tmp = 
			(al0*al0*(3.0 - 2.0*al0) + b1*b1*(b1 - 3.0*al0) 
			+ b2*b2*(b2 - 3.0*al0) + b3*b3*(b3 - 3.0*al0))
			/pr;
	}
      
	double V = 0.0;  
	if (al <= 0.5)
	{
		V = tmp*dr0;
	}
	else
	{
		V = (1.0 - tmp)*dr0;
	}

	return V;
}


double VOF_PLIC::calcV2(lexer *p)
{
	double M1c1 = nx(i,j,k)*p->DXN[IP];
	double m2c2 = ny(i,j,k)*p->DYN[JP];
	double m3c3 = nz(i,j,k)*p->DZN[KP];
	
	double alpha_max = M1c1 + m2c2 + m3c3;
	
	if (fabs(alpha_max) > 1e-6)
	{
		double V = 
			1.0/(6.0*M1c1*m2c2*m3c3)*
			(
				pow(alpha(i,j,k),3.0) 
				- F3(alpha(i,j,k) - M1c1) 
				- F3(alpha(i,j,k) - m2c2) 
				- F3(alpha(i,j,k) - m3c3)
				+ F3(alpha(i,j,k) - alpha_max + M1c1) 
				+ F3(alpha(i,j,k) - alpha_max + m2c2) 
				+ F3(alpha(i,j,k) - alpha_max + m3c3)
			);

		return V;
	}
	else
	{
		return 0.0;
	}
}


double VOF_PLIC::F3(const double& x)
{	
	return (x <= 0.0) ? 0.0 : pow(x,3.0);
}

double VOF_PLIC::Volume_flow(lexer* p, double& nx, double& ny, double& nz, double& alpha, int sweepnum, int dirnum)
{
    double retval, al,alphamax;
    
    double n_x,n_y,n_z,ret;
    n_x=abs(nx);
    n_y=abs(ny);
    n_z=abs(nz);
    
    //sorting by n_i * d_i
    double mx,my,mz,dx,dy,dz,d_x,d_y,d_z;
    d_x=p->DXN[IP];
    d_y=p->DYN[JP];
    d_z=p->DZN[KP];
    
    mx=nx;
    my=ny;
    mz=nz;
    /*
    if(n_x*d_x<=n_y*d_y)
    {
        mx=n_x;
        dx=d_x;
        my=n_y;
        dy=d_y;
	}
    else
    {
        mx=n_y;
        dx=d_y;
        my=n_x;
        dy=d_x;
    }
    
    if(n_z*d_z>=my*dy)
    {
        mz=n_z;
        dz=d_z;
    }
    else if(n_z*d_z<mx*dx)
    {
        mz=my;
        dz=dy;
        my=mx;
        dy=dx;
        mx=n_z;
        dx=d_z;
    }
    else
    {
        mz=my;
        dz=dy;
        my=n_z;
        dy=d_z;
    }
    */
    if (dirnum==2)
    {
        al = alpha;
        
        alphamax=mx*p->DXN[IP1]+my*p->DYN[JP]+mz*p->DZN[KP];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IP1])*(al-mx*p->DXN[IP1])*(al-mx*p->DXN[IP1])*(al-mx*p->DXN[IP1])
                -Heavistep(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])
                -Heavistep(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])
                +Heavistep(al-alphamax+mx*p->DXN[IP1])*(al-alphamax+mx*p->DXN[IP1])*(al-alphamax+mx*p->DXN[IP1])*(al-alphamax+mx*p->DXN[IP1])
                +Heavistep(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])
                +Heavistep(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])
                );
    }
    if (sweepnum==0 && dirnum==3)
    {
        al = alpha-nx*p->DXN[IP];
        
        alphamax=mx*p->DXN[IP1]+my*p->DYN[JP]+mz*p->DZN[KP];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IP1])*(al-mx*p->DXN[IP1])*(al-mx*p->DXN[IP1])*(al-mx*p->DXN[IP1])
                -Heavistep(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])
                -Heavistep(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])
                +Heavistep(al-alphamax+mx*p->DXN[IP1])*(al-alphamax+mx*p->DXN[IP1])*(al-alphamax+mx*p->DXN[IP1])*(al-alphamax+mx*p->DXN[IP1])
                +Heavistep(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])
                +Heavistep(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])
                );
    }
    if (sweepnum==0 && dirnum==1)
    {
        al = alpha-nx*(-1.0)*p->DXN[IP];
        alphamax=mx*p->DXN[IM1]+my*p->DYN[JP]+mz*p->DZN[KP];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IM1])*(al-mx*p->DXN[IM1])*(al-mx*p->DXN[IM1])*(al-mx*p->DXN[IM1])
                -Heavistep(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])
                -Heavistep(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])
                +Heavistep(al-alphamax+mx*p->DXN[IM1])*(al-alphamax+mx*p->DXN[IM1])*(al-alphamax+mx*p->DXN[IM1])*(al-alphamax+mx*p->DXN[IM1])
                +Heavistep(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])
                +Heavistep(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])
                );
    }
    if (sweepnum==1 && dirnum==3)
    {
        al = alpha-ny*p->DYN[JP];
        alphamax=mx*p->DXN[IP]+my*p->DYN[JP1]+mz*p->DZN[KP];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])
                -Heavistep(al-my*p->DYN[JP1])*(al-my*p->DYN[JP1])*(al-my*p->DYN[JP1])*(al-my*p->DYN[JP1])
                -Heavistep(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])
                +Heavistep(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])
                +Heavistep(al-alphamax+my*p->DYN[JP1])*(al-alphamax+my*p->DYN[JP1])*(al-alphamax+my*p->DYN[JP1])*(al-alphamax+my*p->DYN[JP1])
                +Heavistep(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])
                );
    }
    if (sweepnum==1 && dirnum==1)
    {
        al = alpha-ny*(-1.0)*p->DYN[JP];
        alphamax=mx*p->DXN[IP]+my*p->DYN[JM1]+mz*p->DZN[KP];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])
                -Heavistep(al-my*p->DYN[JM1])*(al-my*p->DYN[JM1])*(al-my*p->DYN[JM1])*(al-my*p->DYN[JM1])
                -Heavistep(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])*(al-mz*p->DZN[KP])
                +Heavistep(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])
                +Heavistep(al-alphamax+my*p->DYN[JM1])*(al-alphamax+my*p->DYN[JM1])*(al-alphamax+my*p->DYN[JM1])*(al-alphamax+my*p->DYN[JM1])
                +Heavistep(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])*(al-alphamax+mz*p->DZN[KP])
                );
    }
    if (sweepnum==2 && dirnum==3)
    {
        al = alpha-nz*p->DZN[KP];
        alphamax=mx*p->DXN[IP]+my*p->DYN[JP]+mz*p->DZN[KP1];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])
                -Heavistep(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])
                -Heavistep(al-mz*p->DZN[KP1])*(al-mz*p->DZN[KP1])*(al-mz*p->DZN[KP1])*(al-mz*p->DZN[KP1])
                +Heavistep(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])
                +Heavistep(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])
                +Heavistep(al-alphamax+mz*p->DZN[KP1])*(al-alphamax+mz*p->DZN[KP1])*(al-alphamax+mz*p->DZN[KP1])*(al-alphamax+mz*p->DZN[KP1])
                );
    }
    if (sweepnum==2 && dirnum==1)
    {
        al = alpha-nz*(-1.0)*p->DZN[KP];
        alphamax=mx*p->DXN[IP]+my*p->DYN[JP]+mz*p->DZN[KM1];
        retval = 1/(6*mx*my*mz)
                *(al*al*al
                -Heavistep(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])*(al-mx*p->DXN[IP])
                -Heavistep(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])*(al-my*p->DYN[JP])
                -Heavistep(al-mz*p->DZN[KM1])*(al-mz*p->DZN[KM1])*(al-mz*p->DZN[KM1])*(al-mz*p->DZN[KM1])
                +Heavistep(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])*(al-alphamax+mx*p->DXN[IP])
                +Heavistep(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])*(al-alphamax+my*p->DYN[JP])
                +Heavistep(al-alphamax+mz*p->DZN[KM1])*(al-alphamax+mz*p->DZN[KM1])*(al-alphamax+mz*p->DZN[KM1])*(al-alphamax+mz*p->DZN[KM1])
                );
    }
    if(sqrt(mx*mx+my*my+mz*mz)<0.9 || sqrt(mx*mx+my*my+mz*mz)>1.1)
        cout<<"Nonuniform normvec in Vol FLux!"<<endl;
    //if(al<0.0)
        //cout<<"Neg al problem in Vol Flux!"<<endl;
    if(retval<0.0)
    {
        //cout<<"Neg Vol Flux Alarm!: "<<"sweepnum: "<<sweepnum<<" dirnum: "<<dirnum<<" amount per cellvolume "<<retval/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP])<<endl;
        retval=0.0;
    }
    if(retval<=0.0)
        retval=0.0;
    return retval;
}

double VOF_PLIC::Heavistep(double x)
{
    double retval;
    if(x>0.0)
        retval = 1.0;
    else
        retval = 0.0;
    return retval;
} 