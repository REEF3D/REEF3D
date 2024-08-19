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
Author: Fabian Knoblauch
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

void VOF_PLIC::advectPlane_alt
(
	fdm* a, 
	lexer* p, 
	int sweep
)
{
	double ds_p, ds_m, ds_c;
    double nx_mod, ny_mod, nz_mod,r0_mod;
    double r0 = alpha(i,j,k);
    if(sweep==0)
    {
        ds_p=a->u(i,j,k)*p->dt;
        ds_m=a->u(i-1,j,k)*p->dt;
        ny_mod=ny(i,j,k);
        nz_mod=nz(i,j,k);
        nx_mod=nx(i,j,k)*p->DXN[IP]/(p->DXN[IP]+ds_p-ds_m);
        ds_c=(a->u(i,j,k)+a->u(i-1,j,k))/2.0*p->dt;
        r0_mod=nx_mod*nx_mod*r0+ny_mod*ny_mod*r0+nz_mod*nz_mod*r0+nx_mod*ds_c;
        if(ds_p>0.0)
        {
            double r0_p,V_w_p;
            r0_p=nx_mod*(p->DXN[IP]/2.0+ds_p/2.0)-r0_mod;
            V_w_p=calculateVolume(nx_mod,ny_mod,nz_mod,ds_p,p->DYN[JP],p->DZN[KP],r0_p); //outputs Volume fraction [0;1] of cellsize
            V_w_p=V_w_p*ds_p*p->DYN[JP]*p->DZN[KP]; //translation into absolute Volume
            V_w_update(i+1,j,k)+=V_w_p;
            V_w_update(i,j,k)-=V_w_p;
            V_a_update(i+1,j,k)+=(ds_p*p->DYN[JP]*p->DZN[KP])-V_w_p;
            V_a_update(i,j,k)-=(ds_p*p->DYN[JP]*p->DZN[KP])-V_w_p;
            
            if(vof_old(i+1,j,k)>0.9999)
            {
                V_w_update(i+1,j,k)-=ds_p*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i+1,j,k)<0.0001)
            {
                V_a_update(i+1,j,k)-=ds_p*p->DYN[JP]*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i+1,j,k)>0.9999)
            {
                V_w_update(i,j,k)+=fabs(ds_p)*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i+1,j,k)<0.0001)
            {
                V_a_update(i,j,k)+=fabs(ds_p)*p->DYN[JP]*p->DZN[KP];
            }
        }
        
        if(ds_m<0.0)
        {
            double r0_m,V_w_m;
            r0_m=nx_mod*(-p->DXN[IP]/2.0+ds_m/2.0)-r0_mod;
            V_w_m=calculateVolume(nx_mod,ny_mod,nz_mod,fabs(ds_m),p->DYN[JP],p->DZN[KP],r0_m); 
            V_w_m=V_w_m*fabs(ds_m)*p->DYN[JP]*p->DZN[KP]; //
            V_w_update(i-1,j,k)+=V_w_m;
            V_w_update(i,j,k)-=V_w_m;
            V_a_update(i-1,j,k)+=(fabs(ds_m)*p->DYN[JP]*p->DZN[KP])-V_w_m;
            V_a_update(i,j,k)-=(fabs(ds_m)*p->DYN[JP]*p->DZN[KP])-V_w_m;
            
            if(vof_old(i-1,j,k)>0.9999)
            {
                V_w_update(i-1,j,k)-=fabs(ds_m)*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i-1,j,k)<0.0001)
            {
                V_a_update(i-1,j,k)-=fabs(ds_m)*p->DYN[JP]*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i-1,j,k)>0.9999)
            {
                V_w_update(i,j,k)+=ds_m*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i-1,j,k)<0.0001)
            {
                V_a_update(i,j,k)+=ds_m*p->DYN[JP]*p->DZN[KP];
            }
        }
        
    }
    else if(sweep==1)
    {
        ds_p=a->v(i,j,k)*p->dt;
        ds_m=a->v(i,j-1,k)*p->dt;
        ny_mod=ny(i,j,k)*p->DYN[JP]/(p->DYN[JP]+ds_p-ds_m);
        nz_mod=nz(i,j,k);
        nx_mod=nx(i,j,k);
        ds_c=(a->v(i,j,k)+a->v(i,j-1,k))/2.0*p->dt;
        r0_mod=nx_mod*nx_mod*r0+ny_mod*ny_mod*r0+nz_mod*nz_mod*r0+ny_mod*ds_c;
        if(ds_p>0.0)
        {
            double r0_p,V_w_p;
            r0_p=ny_mod*(p->DYN[JP]/2.0+ds_p/2.0)-r0_mod;
            V_w_p=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],ds_p,p->DZN[KP],r0_p); //outputs Volume fraction [0;1] of cellsize
            V_w_p=V_w_p*p->DXN[IP]*ds_p*p->DZN[KP]; //translation into absolute Volume
            V_w_update(i,j+1,k)+=V_w_p;
            V_w_update(i,j,k)-=V_w_p;
            V_a_update(i,j+1,k)+=(p->DXN[IP]*ds_p*p->DZN[KP])-V_w_p;
            V_a_update(i,j,k)-=(p->DXN[IP]*ds_p*p->DZN[KP])-V_w_p;
            
            if(vof_old(i,j+1,k)>0.9999)
            {
                V_w_update(i,j+1,k)-=p->DXN[IP]*ds_p*p->DZN[KP];
            }
            else if(vof_old(i,j+1,k)<0.0001)
            {
                V_a_update(i,j+1,k)-=p->DXN[IP]*ds_p*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i,j+1,k)>0.9999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*fabs(ds_p)*p->DZN[KP];
            }
            else if(vof_old(i,j+1,k)<0.0001)
            {
                V_a_update(i,j,k)+=p->DXN[IP]*fabs(ds_p)*p->DZN[KP];
            }
        }
        
        if(ds_m<0.0)
        {
            double r0_m,V_w_m;
            r0_m=ny_mod*(-p->DYN[JP]/2.0+ds_m/2.0)-r0_mod;
            V_w_m=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],fabs(ds_m),p->DZN[KP],r0_m); 
            V_w_m=V_w_m*p->DXN[IP]*fabs(ds_m)*p->DZN[KP]; 
            V_w_update(i,j-1,k)+=V_w_m;
            V_w_update(i,j,k)-=V_w_m;
            V_a_update(i,j-1,k)+=(p->DXN[IP]*fabs(ds_m)*p->DZN[KP])-V_w_m;
            V_a_update(i,j,k)-=(p->DXN[IP]*fabs(ds_m)*p->DZN[KP])-V_w_m;
            
            if(vof_old(i,j-1,k)>0.9999)
            {
                V_w_update(i,j-1,k)-=p->DXN[IP]*fabs(ds_m)*p->DZN[KP];
            }
            else if(vof_old(i,j-1,k)<0.0001)
            {
                V_a_update(i,j-1,k)-=p->DXN[IP]*fabs(ds_m)*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i,j-1,k)>0.9999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*ds_m*p->DZN[KP];
            }
            else if(vof_old(i,j-1,k)<0.0001)
            {
                V_a_update(i,j,k)+=p->DXN[IP]*ds_m*p->DZN[KP];
            }
        }
        
    }
    
    else if(sweep==2)
    {
        ds_p=a->w(i,j,k)*p->dt;
        ds_m=a->w(i,j,k-1)*p->dt;
        ny_mod=ny(i,j,k);
        nz_mod=nz(i,j,k)*p->DZN[KP]/(p->DZN[KP]+ds_p-ds_m);
        nx_mod=nx(i,j,k);
        ds_c=(a->w(i,j,k)+a->w(i,j,k-1))/2.0*p->dt;
        r0_mod=nx_mod*nx_mod*r0+ny_mod*ny_mod*r0+nz_mod*nz_mod*r0+nz_mod*ds_c;
        if(ds_p>0.0)
        {
            double r0_p,V_w_p;
            r0_p=nz_mod*(p->DZN[KP]/2.0+ds_p/2.0)-r0_mod;
            V_w_p=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],p->DYN[JP],ds_p,r0_p); //outputs Volume fraction [0;1] of cellsize
            V_w_p=V_w_p*p->DXN[IP]*p->DYN[JP]*ds_p; //translation into absolute Volume
            V_w_update(i,j,k+1)+=V_w_p;
            V_w_update(i,j,k)-=V_w_p;
            V_a_update(i,j,k+1)+=(p->DXN[IP]*p->DYN[JP]*ds_p)-V_w_p;
            V_a_update(i,j,k)-=(p->DXN[IP]*p->DYN[JP]*ds_p)-V_w_p;
            
            if(vof_old(i,j,k+1)>0.9999)
            {
                V_w_update(i,j,k+1)-=p->DXN[IP]*p->DYN[JP]*ds_p;
            }
            else if(vof_old(i,j,k+1)<0.0001)
            {
                V_a_update(i,j,k+1)-=p->DXN[IP]*p->DYN[JP]*ds_p;
            }
        }
        else
        {
            if(vof_old(i,j,k+1)>0.9999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*fabs(ds_p);
            }
            else if(vof_old(i,j,k+1)<0.0001)
            {
                V_a_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*fabs(ds_p);
            }
        }
        
        if(ds_m<0.0)
        {
            double r0_m,V_w_m;
            r0_m=nz_mod*(-p->DZN[KP]/2.0+ds_m/2.0)-r0_mod;
            V_w_m=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],p->DYN[JP],fabs(ds_m),r0_m); 
            V_w_m=V_w_m*p->DXN[IP]*p->DYN[JP]*fabs(ds_m); 
            V_w_update(i,j,k-1)+=V_w_m;
            V_w_update(i,j,k)-=V_w_m;
            V_a_update(i,j,k-1)+=(p->DXN[IP]*p->DYN[JP]*fabs(ds_m))-V_w_m;
            V_a_update(i,j,k)-=(p->DXN[IP]*p->DYN[JP]*fabs(ds_m))-V_w_m;
            
            if(vof_old(i,j,k-1)>0.9999)
            {
                V_w_update(i,j,k-1)-=p->DXN[IP]*p->DYN[JP]*fabs(ds_m);
            }
            else if(vof_old(i,j,k-1)<0.0001)
            {
                V_a_update(i,j,k-1)-=p->DXN[IP]*p->DYN[JP]*fabs(ds_m);
            }
        }
        else
        {
            if(vof_old(i,j,k-1)>0.9999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*ds_m;
            }
            else if(vof_old(i,j,k-1)<0.0001)
            {
                V_a_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*ds_m;
            }
        }
    }
}

double VOF_PLIC::calculateVolume(double n_a, double n_b, double n_c, double d_a, double d_b, double d_c, double r0)
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
    
    if(r0>=0)
    {
        V0=(0.5-V)+0.5;
    }
    else
    {
        V0=-(0.5-V)+0.5;
    }
    
    if(V0<0.0-1E6)
        cout<<"neg VO output"<<endl;
    if(V0>1.0+1E6)
        cout<<"too hight V0 output"<<endl;
        
    return V0;
}