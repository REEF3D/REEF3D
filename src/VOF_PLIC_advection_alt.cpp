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
            
            if(vof_old(i+1,j,k)>0.999)
            {
                V_w_update(i+1,j,k)-=ds_p*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i+1,j,k)<0.001)
            {
                V_a_update(i+1,j,k)-=ds_p*p->DYN[JP]*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i+1,j,k)>0.999)
            {
                V_w_update(i,j,k)+=fabs(ds_p)*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i+1,j,k)<0.001)
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
            
            if(vof_old(i-1,j,k)>0.999)
            {
                V_w_update(i-1,j,k)-=fabs(ds_m)*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i-1,j,k)<0.001)
            {
                V_a_update(i-1,j,k)-=fabs(ds_m)*p->DYN[JP]*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i-1,j,k)>0.999)
            {
                V_w_update(i,j,k)+=ds_m*p->DYN[JP]*p->DZN[KP];
            }
            else if(vof_old(i-1,j,k)<0.001)
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
            
            if(vof_old(i,j+1,k)>0.999)
            {
                V_w_update(i,j+1,k)-=p->DXN[IP]*ds_p*p->DZN[KP];
            }
            else if(vof_old(i,j+1,k)<0.001)
            {
                V_a_update(i,j+1,k)-=p->DXN[IP]*ds_p*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i,j+1,k)>0.999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*fabs(ds_p)*p->DZN[KP];
            }
            else if(vof_old(i,j+1,k)<0.001)
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
            
            if(vof_old(i,j-1,k)>0.999)
            {
                V_w_update(i,j-1,k)-=p->DXN[IP]*fabs(ds_m)*p->DZN[KP];
            }
            else if(vof_old(i,j-1,k)<0.001)
            {
                V_a_update(i,j-1,k)-=p->DXN[IP]*fabs(ds_m)*p->DZN[KP];
            }
        }
        else
        {
            if(vof_old(i,j-1,k)>0.999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*ds_m*p->DZN[KP];
            }
            else if(vof_old(i,j-1,k)<0.001)
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
            
            if(vof_old(i,j,k+1)>0.999)
            {
                V_w_update(i,j,k+1)-=p->DXN[IP]*p->DYN[JP]*ds_p;
            }
            else if(vof_old(i,j,k+1)<0.001)
            {
                V_a_update(i,j,k+1)-=p->DXN[IP]*p->DYN[JP]*ds_p;
            }
        }
        else
        {
            if(vof_old(i,j,k+1)>0.999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*fabs(ds_p);
            }
            else if(vof_old(i,j,k+1)<0.001)
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
            
            if(vof_old(i,j,k-1)>0.999)
            {
                V_w_update(i,j,k-1)-=p->DXN[IP]*p->DYN[JP]*fabs(ds_m);
            }
            else if(vof_old(i,j,k-1)<0.001)
            {
                V_a_update(i,j,k-1)-=p->DXN[IP]*p->DYN[JP]*fabs(ds_m);
            }
        }
        else
        {
            if(vof_old(i,j,k-1)>0.999)
            {
                V_w_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*ds_m;
            }
            else if(vof_old(i,j,k-1)<0.001)
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
    
    if(V0<0.0)
        cout<<"neg VO output"<<endl;
    if(V0>1.0)
        cout<<"too hight V0 output"<<endl;
        
    return V0;
}


void VOF_PLIC::advectPlane_altFlux
(
	fdm* a, 
	lexer* p, 
    double Q1,
    double Q2,
	int sweep
)
{
	double ds_p, ds_m, ds_c;
    double nx_mod, ny_mod, nz_mod,r0_mod, nmodsum;
    double r0 = alpha(i,j,k);
    if(sweep==0)
    {
        ds_p=Q2*p->DXN[IP];
        ds_m=Q1*p->DXN[IP];
        ny_mod=ny(i,j,k);
        nz_mod=nz(i,j,k);
        nx_mod=nx(i,j,k)*p->DXN[IP]/(p->DXN[IP]+ds_p-ds_m);
        nmodsum=sqrt(nx_mod*nx_mod+ny_mod*ny_mod+nz_mod*nz_mod);
        nx_mod=nx_mod/nmodsum;
        ny_mod=ny_mod/nmodsum;
        nz_mod=nz_mod/nmodsum;
        r0=r0/nmodsum;
        if(nx_mod>=0.0)
            ds_c=ds_m;
        else
            ds_c=ds_p;
        r0_mod=nx_mod*nx_mod*r0+ny_mod*ny_mod*r0+nz_mod*nz_mod*r0+nx_mod*ds_c;
        if(ds_p>1E-10)
        {
            double r0_p,V_w_p;
            r0_p=nx_mod*(p->DXN[IP]/2.0+ds_p/2.0)-r0_mod;
            V_w_p=calculateVolume(nx_mod,ny_mod,nz_mod,ds_p,p->DYN[JP],p->DZN[KP],r0_p); //outputs Volume fraction [0;1] of cellsize
            V_w_p=V_w_p*ds_p*p->DYN[JP]*p->DZN[KP]; //translation into absolute Volume
            V_w_update(i+1,j,k)+=V_w_p;
            V_w_update(i,j,k)-=V_w_p;
        }
        
        if(ds_m<-1E-10)
        {
            double r0_m,V_w_m;
            r0_m=nx_mod*(-p->DXN[IP]/2.0+ds_m/2.0)-r0_mod;
            V_w_m=calculateVolume(nx_mod,ny_mod,nz_mod,fabs(ds_m),p->DYN[JP],p->DZN[KP],r0_m); 
            V_w_m=V_w_m*fabs(ds_m)*p->DYN[JP]*p->DZN[KP]; //
            V_w_update(i-1,j,k)+=V_w_m;
            V_w_update(i,j,k)-=V_w_m;
        }
        
        nx(i,j,k)=nx_mod;
        ny(i,j,k)=ny_mod;
        nz(i,j,k)=nz_mod;
        alpha(i,j,k)=r0_mod;
        
    }
    else if(sweep==1)
    {
        ds_p=Q2*p->DYN[JP];
        ds_m=Q1*p->DYN[JP];
        ny_mod=ny(i,j,k)*p->DYN[JP]/(p->DYN[JP]+ds_p-ds_m);
        nz_mod=nz(i,j,k);
        nx_mod=nx(i,j,k);
        nmodsum=sqrt(nx_mod*nx_mod+ny_mod*ny_mod+nz_mod*nz_mod);
        nx_mod=nx_mod/nmodsum;
        ny_mod=ny_mod/nmodsum;
        nz_mod=nz_mod/nmodsum;
        r0=r0/nmodsum;
        if(ny_mod>=0.0)
            ds_c=ds_m;
        else
            ds_c=ds_p;
        r0_mod=nx_mod*nx_mod*r0+ny_mod*ny_mod*r0+nz_mod*nz_mod*r0+ny_mod*ds_c;
        if(ds_p>1E-10)
        {
            double r0_p,V_w_p;
            r0_p=ny_mod*(p->DYN[JP]/2.0+ds_p/2.0)-r0_mod;
            V_w_p=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],ds_p,p->DZN[KP],r0_p); //outputs Volume fraction [0;1] of cellsize
            V_w_p=V_w_p*p->DXN[IP]*ds_p*p->DZN[KP]; //translation into absolute Volume
            V_w_update(i,j+1,k)+=V_w_p;
            V_w_update(i,j,k)-=V_w_p;
        }   
        
        if(ds_m<-1E-10)
        {
            double r0_m,V_w_m;
            r0_m=ny_mod*(-p->DYN[JP]/2.0+ds_m/2.0)-r0_mod;
            V_w_m=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],fabs(ds_m),p->DZN[KP],r0_m); 
            V_w_m=V_w_m*p->DXN[IP]*fabs(ds_m)*p->DZN[KP]; 
            V_w_update(i,j-1,k)+=V_w_m;
            V_w_update(i,j,k)-=V_w_m;
        }
        nx(i,j,k)=nx_mod;
        ny(i,j,k)=ny_mod;
        nz(i,j,k)=nz_mod;
        alpha(i,j,k)=r0_mod;
    }
    
    else if(sweep==2)
    {
        ds_p=Q2*p->DZN[KP];
        ds_m=Q1*p->DZN[KP];
        ny_mod=ny(i,j,k);
        nz_mod=nz(i,j,k)*p->DZN[KP]/(p->DZN[KP]+ds_p-ds_m);
        nx_mod=nx(i,j,k);
        nmodsum=sqrt(nx_mod*nx_mod+ny_mod*ny_mod+nz_mod*nz_mod);
        nx_mod=nx_mod/nmodsum;
        ny_mod=ny_mod/nmodsum;
        nz_mod=nz_mod/nmodsum;
        r0=r0/nmodsum;
        if(nz_mod>=0.0)
            ds_c=ds_m;
        else
            ds_c=ds_p;
        r0_mod=nx_mod*nx_mod*r0+ny_mod*ny_mod*r0+nz_mod*nz_mod*r0+nz_mod*ds_c;
        if(ds_p>0.0)
        {
            double r0_p,V_w_p;
            r0_p=nz_mod*(p->DZN[KP]/2.0+ds_p/2.0)-r0_mod;
            V_w_p=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],p->DYN[JP],ds_p,r0_p); //outputs Volume fraction [0;1] of cellsize
            V_w_p=V_w_p*p->DXN[IP]*p->DYN[JP]*ds_p; //translation into absolute Volume
            V_w_update(i,j,k+1)+=V_w_p;
            V_w_update(i,j,k)-=V_w_p;
        }
        if(ds_m<0.0)
        {
            double r0_m,V_w_m;
            r0_m=nz_mod*(-p->DZN[KP]/2.0+ds_m/2.0)-r0_mod;
            V_w_m=calculateVolume(nx_mod,ny_mod,nz_mod,p->DXN[IP],p->DYN[JP],fabs(ds_m),r0_m); 
            V_w_m=V_w_m*p->DXN[IP]*p->DYN[JP]*fabs(ds_m); 
            V_w_update(i,j,k-1)+=V_w_m;
            V_w_update(i,j,k)-=V_w_m;
        }
        nx(i,j,k)=nx_mod;
        ny(i,j,k)=ny_mod;
        nz(i,j,k)=nz_mod;
        alpha(i,j,k)=r0_mod;
    }
}

void VOF_PLIC::advectPlane_sweepless
(
	fdm* a, 
	lexer* p
)
{
    double velx,vely,velz,dsx,dsy,dsz,r0mod;
    
    if(nx(i,j,k)>0.0)
    {
        velx=a->u(i-1,j,k);
    }
    else if(nx(i,j,k)<0.0)
    {
        velx=a->u(i,j,k);
    }
    else
    {
        velx=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    }
    dsx=velx*p->dt;
    if(fabs(dsx)>0.5*p->DXN[IP])
            cout<<"dsx faster than CFL"<<endl;
    
    if(ny(i,j,k)>0.0)
    {
        vely=a->v(i,j-1,k);
    }
    else if(ny(i,j,k)<0.0)
    {
        vely=a->v(i,j,k);
    }   
    else
    {
        vely=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    }
    dsy=vely*p->dt;
        if(fabs(dsy)>0.5*p->DYN[JP])
            cout<<"dy faster than CFL"<<endl;
    
    if(nz(i,j,k)>0.0)
    {
        velz=a->w(i,j,k-1);
    }   
    else if(nz(i,j,k)<0.0)
    {
        velz=a->w(i,j,k);
    }   
    else
    {
        velz=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    }
    dsz=velz*p->dt;
         if(fabs(dsz)>0.5*p->DZN[KP])
            cout<<"dsz faster than CFL"<<endl;
    
    r0mod=alpha(i,j,k)*(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k))+nx(i,j,k)*dsx+ny(i,j,k)*dsy+nz(i,j,k)*dsz;
    //cout<<"alpha:"<<alpha(i,j,k)<<" r0_mod: "<<r0mod<<endl;

    // 2D only
    if (fabs(velx)>=1E-6)
    {   
        if(dsx>0.0)
        {
            if(a->vof(i+1,j,k)<0.999)
            {
                double r0xdiff, transfVolx, recheck, scaledVol, absVol;
                r0xdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0xdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0xdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                    V_w_update(i+1,j,k)+=absVol;
                    V_w_update(i,j,k)-=absVol;
                }
                if(a->vof(i-1,j,k)>0.999)
                    V_w_update(i,j,k)+=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
            }
            else
            {
                double r0xdiff, transfVolx, recheck, scaledVol, absVol, takenVol;
                r0xdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0xdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0xdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                    takenVol=(1.0-scaledVol)*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                    V_w_update(i+1,j,k)-=takenVol;
                    V_w_update(i,j,k)-=absVol;
                }
                else
                    V_w_update(i,j,k)-=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
            }
        }
        else
        {
            if(a->vof(i-1,j,k)<0.999)
            {
                double r0xdiff, transfVolx, recheck, scaledVol, absVol;
                r0xdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0xdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0xdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                    V_w_update(i-1,j,k)+=absVol;
                    V_w_update(i,j,k)-=absVol;
                }
                if(a->vof(i+1,j,k)>0.999)
                    V_w_update(i,j,k)+=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
            }
            else
            {
                double r0xdiff, transfVolx, recheck, scaledVol, absVol, takenVol;
                r0xdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0xdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0xdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                    takenVol=(1.0-scaledVol)*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                    V_w_update(i-1,j,k)-=takenVol;
                    V_w_update(i,j,k)-=absVol;
                }
                else
                    V_w_update(i,j,k)-=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
            }
        }
    }
    
    if (fabs(velz)>=1E-6)
    {   
        if(dsz>0.0)
        {
            if(a->vof(i,j,k+1)<0.999)
            {
                double r0zdiff, transfVolz, recheck, scaledVol, absVol;
                r0zdiff=nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*p->DXN[IP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0zdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz), r0zdiff);
                    absVol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                    V_w_update(i,j,k-1)+=absVol;
                    V_w_update(i,j,k)-=absVol;
                }
                if(a->vof(i,j,k-1)>0.999)
                    V_w_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
            }
            else
            {
                double r0zdiff, transfVolz, recheck, scaledVol, absVol, takenVol;
                r0zdiff=nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*p->DXN[IP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0zdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0zdiff);
                    absVol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                    takenVol=(1.0-scaledVol)*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                    V_w_update(i,j,k+1)-=takenVol;
                    V_w_update(i,j,k)-=absVol;
                }
                else
                    V_w_update(i,j,k)-=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
            }
        }
        else
        {
            if(a->vof(i,j,k-1)<0.999)
            {
                double r0zdiff, transfVolz, recheck, scaledVol, absVol;
                r0zdiff=nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*p->DXN[IP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0zdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz), r0zdiff);
                    absVol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                    V_w_update(i,j,k-1)+=absVol;
                    V_w_update(i,j,k)-=absVol;
                }
                if(a->vof(i,j,k+1)>0.999)
                    V_w_update(i,j,k)+=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
            }
            else
            {
                double r0zdiff, transfVolz, recheck, scaledVol, absVol, takenVol;
                r0zdiff=nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*p->DXN[IP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0zdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0zdiff);
                    absVol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                    takenVol=(1.0-scaledVol)*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                    V_w_update(i,j,k-1)-=takenVol;
                    V_w_update(i,j,k)-=absVol;
                }
                else
                    V_w_update(i,j,k)-=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
            }
        }
        
    }
    
    if (fabs(velz)>=1E-6 && fabs(velx)>=1E-6)
    {
        if(dsx>0.0 && dsz>0.0)
        {
            if(a->vof(i+1,j,k+1)<0.999)
            {
                double r0xzdiff, recheck, scaledVol, absVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i+1,j,k+1)+=absVol;
                    V_w_old(i,j,k)=p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                }
            }
            else
            {
                double r0xzdiff, recheck, scaledVol, absVol, takenVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    takenVol=(1.0-scaledVol)*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i+1,j,k+1)-=takenVol;
                    V_w_old(i,j,k)=0.0;
                }
            }
        }
        else if(dsx>0.0 && dsz<0.0)
        {
            if(a->vof(i+1,j,k-1)<0.999)
            {
                double r0xzdiff, recheck, scaledVol, absVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i+1,j,k-1)+=absVol;
                    V_w_old(i,j,k)=p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                }
            }
            else
            {
                double r0xzdiff, recheck, scaledVol, absVol, takenVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    takenVol=(1.0-scaledVol)*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i+1,j,k-1)-=takenVol;
                    V_w_old(i,j,k)=0.0;
                }
            }
        }
        else if(dsx<0.0 && dsz>0.0)
        {
            if(a->vof(i-1,j,k+1)<0.999)
            {
                double r0xzdiff, recheck, scaledVol, absVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i-1,j,k+1)+=absVol;
                    V_w_old(i,j,k)=p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                }
            }
            else
            {
                double r0xzdiff, recheck, scaledVol, absVol, takenVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    takenVol=(1.0-scaledVol)*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i-1,j,k+1)-=takenVol;
                    V_w_old(i,j,k)=0.0;
                }
            }
        }
        else
        {
            if(a->vof(i-1,j,k-1)<0.999)
            {
                double r0xzdiff, recheck, scaledVol, absVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    absVol=scaledVol*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i-1,j,k-1)+=absVol;
                    V_w_old(i,j,k)=p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                }
            }
            else
            {
                double r0xzdiff, recheck, scaledVol, absVol, takenVol;
                r0xzdiff=nx(i,j,k)*copysign(0.5*p->DXN[IP]+0.5*fabs(dsx),dsx)+nz(i,j,k)*copysign(0.5*p->DZN[KP]+0.5*fabs(dsz),dsz)-r0mod;
                recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0xzdiff);
                if(recheck>0.0)
                {
                    scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],fabs(dsz),r0xzdiff);
                    takenVol=(1.0-scaledVol)*fabs(dsx)*p->DYN[JP]*fabs(dsz);
                    V_w_update(i-1,j,k-1)-=takenVol;
                    V_w_old(i,j,k)=0.0;
                }
            }
        }
        
    }
    
}

void VOF_PLIC::advectWater_sweepless
(
	fdm* a, 
	lexer* p
)
{
    double velx,vely,velz,dsx,dsy,dsz;
    
    if(a->vof(i+1,j,k)>a->vof(i-1,j,k))
    {
        velx=a->u(i,j,k);
    }   
    else if(a->vof(i+1,j,k)<a->vof(i-1,j,k))
    {
        velx=a->u(i-1,j,k);
    }   
    else
    {
        velx=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    }
    
    if(a->vof(i,j+1,k)>a->vof(i,j-1,k))
    {
        vely=a->v(i,j,k);
    }
    else if(a->vof(i,j+1,k)<a->vof(i,j-1,k))
    {
        vely=a->v(i-1,j,k);
    }
    else
    {
        vely=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    }
    
    if(a->vof(i,j,k+1)>a->vof(i,j,k-1))
    {
        velz=a->w(i,j,k);
    }
    else if(a->vof(i,j,k+1)<a->vof(i,j,k-1))
    {
        velz=a->w(i,j,k-1);
    }    
    else
    {
        vely=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    }
    
    
    dsx=velx*p->dt;
    if(fabs(dsx)>0.5*p->DXN[IP])
            cout<<"dsx faster than CFL"<<endl;
    dsy=vely*p->dt;
    if(fabs(dsy)>0.5*p->DYN[JP])
            cout<<"dsy faster than CFL"<<endl;
    dsz=velz*p->dt;
    if(fabs(dsz)>0.5*p->DZN[KP])
            cout<<"dsz faster than CFL"<<endl;
            
    // 2D
    if (fabs(velx)>=1E-6)
    {   
        if(dsx>0.0)
        {
            V_w_update(i+1,j,k)+=fabs(dsx)*p->DYN[IP]*p->DZN[KP];
            V_w_update(i,j,k)-=fabs(dsx)*p->DYN[IP]*p->DZN[KP];
        }
        else
        {
            V_w_update(i-1,j,k)+=fabs(dsx)*p->DYN[IP]*p->DZN[KP];
            V_w_update(i,j,k)-=fabs(dsx)*p->DYN[IP]*p->DZN[KP];
        }
    }
    if (fabs(velz)>=1E-6)
    {   
        if(dsz>0.0)
        {
            V_w_update(i,j,k+1)+=p->DXN[IP]*p->DYN[IP]*fabs(dsz);
            V_w_update(i,j,k)-=p->DXN[IP]*p->DYN[IP]*fabs(dsz);
        }
        else
        {
            V_w_update(i,j,k-1)+=p->DXN[IP]*p->DYN[IP]*fabs(dsz);
            V_w_update(i,j,k)-=p->DXN[IP]*p->DYN[IP]*fabs(dsz);
        }
    }
    if (fabs(velx)>=1E-6 && fabs(velz)>=1E-6)
    {
        if(dsx>0.0 && dsz>0.0)
        {
            V_w_update(i+1,j,k+1)+=fabs(dsx)*p->DYN[JP]*fabs(dsz);
            V_w_update(i,j,k)-=fabs(dsx)*p->DYN[JP]*fabs(dsz);
        }
        else if(dsx>0.0 && dsz<0.0)
        {
            V_w_update(i+1,j,k-1)+=fabs(dsx)*p->DYN[JP]*fabs(dsz);
            V_w_update(i,j,k)-=fabs(dsx)*p->DYN[JP]*fabs(dsz);
        }
        else if(dsx<0.0 && dsz>0.0)
        {
            V_w_update(i-1,j,k+1)+=fabs(dsx)*p->DYN[JP]*fabs(dsz);
            V_w_update(i,j,k)-=fabs(dsx)*p->DYN[JP]*fabs(dsz);
        }
        else
        {
            V_w_update(i-1,j,k-1)+=fabs(dsx)*p->DYN[JP]*fabs(dsz);
            V_w_update(i,j,k)-=fabs(dsx)*p->DYN[JP]*fabs(dsz);
        }    
    }
}
    
void VOF_PLIC::advectPlane_Weymouth
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    // 2D
    if(sweep==0)
    {
        if(a->u(i,j,k)>1E-20)
        {
            double dsx, scaledVol, Vol, r0x, recheck;
            dsx=a->u(i,j,k)*p->dt;
            r0x=nx(i,j,k)*(0.5*p->DXN[IP]-0.5*dsx)-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_update(i+1,j,k)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
            else if(nx(i,j,k)<0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_update(i+1,j,k)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
            
        }
        if(a->u(i-1,j,k)<-1E-20)
        {
            double dsx, scaledVol, Vol, r0x,recheck;
            dsx=a->u(i-1,j,k)*p->dt;
            r0x=nx(i,j,k)*(-0.5*p->DXN[IP]+0.5*fabs(dsx))-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_update(i-1,j,k)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
            else if(nx(i,j,k)>0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_update(i-1,j,k)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
        }
    }
    else if(sweep==2)
    {
        if(a->w(i,j,k)>0.0)
        {
            double dsz,scaledVol, Vol, r0z, recheck;
            dsz=a->w(i,j,k)*p->dt;
            r0z=nz(i,j,k)*(0.5*p->DZN[KP]-0.5*dsz)-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*p->DXN[KP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_update(i,j,k+1)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
            else if(nz(i,j,k)<0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_update(i,j,k+1)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double dsz,scaledVol, Vol, r0z,recheck;
            dsz=a->w(i,j,k-1)*p->dt;
            r0z=nz(i,j,k)*(-0.5*p->DZN[KP]+0.5*fabs(dsz))-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*p->DXN[KP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_update(i,j,k-1)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
            else if(nz(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_update(i,j,k-1)+=Vol;
                V_w_update(i,j,k)-=Vol;
            }
            
        }   
    }
}

void VOF_PLIC::advectWater_Weymouth
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    // 2D
    if(sweep==0)
    {
        if(a->u(i,j,k)>1E-20)
        {
            double Vol;
            Vol=a->u(i,j,k)*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_update(i+1,j,k)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
        if(a->u(i-1,j,k)<-1E-20)
        {
            double Vol;
            Vol=fabs(a->u(i-1,j,k))*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_update(i-1,j,k)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
    }
    else if(sweep==2)
    {
        if(a->w(i,j,k)>1E20)
        {
            double Vol;
            Vol=a->w(i,j,k)*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_update(i,j,k+1)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
        if(a->w(i,j,k-1)<-1E-20)
        {
            double Vol;
            Vol=fabs(a->w(i,j,k-1))*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_update(i,j,k-1)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
    }
}

void VOF_PLIC::advectWater_WeymouthNoS
(
    fdm* a,
    lexer* p
)
{
    // 2D
    
        if(a->u(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->u(i,j,k)*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_update(i+1,j,k)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double Vol;
            Vol=fabs(a->u(i-1,j,k))*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_update(i-1,j,k)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }

   
        if(a->w(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->w(i,j,k)*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_update(i,j,k+1)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double Vol;
            Vol=fabs(a->w(i,j,k-1))*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_update(i,j,k-1)+=Vol;
            V_w_update(i,j,k)-=Vol;
        }
    
}

void VOF_PLIC::advectPlane_MACHO2D
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    // 2D
    if(sweep==0)
    {
        if(a->u(i,j,k)>0.0)
        {
            double dsx, scaledVol, Vol, r0x, recheck;
            dsx=a->u(i,j,k)*p->dt;
            r0x=nx(i,j,k)*(0.5*p->DXN[IP]-0.5*dsx)-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
                V_w_m(i+1,j,k)+=Vol;
                
            }
            else if(nx(i,j,k)<0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
                V_w_m(i+1,j,k)+=Vol;
            }
            
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double dsx, scaledVol, Vol, r0x,recheck;
            dsx=a->u(i-1,j,k)*p->dt;
            r0x=nx(i,j,k)*(-0.5*p->DXN[IP]+0.5*fabs(dsx))-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*fabs(dsx)+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
                V_w_p(i-1,j,k)-=Vol;
            }
            else if(nx(i,j,k)>0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
                V_w_p(i-1,j,k)-=Vol;
            }
        }
    }
    else if(sweep==2)
    {
        if(a->w(i,j,k)>0.0)
        {
            double dsz,scaledVol, Vol, r0z, recheck;
            dsz=a->w(i,j,k)*p->dt;
            r0z=nz(i,j,k)*(0.5*p->DZN[KP]-0.5*dsz)-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*p->DXN[KP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_p(i,j,k)+=Vol;
                V_w_m(i,j,k+1)+=Vol;
            }
            else if(nz(i,j,k)<0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_p(i,j,k)+=Vol;
                V_w_m(i,j,k+1)+=Vol;
            }
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double dsz,scaledVol, Vol, r0z,recheck;
            dsz=a->w(i,j,k-1)*p->dt;
            r0z=nz(i,j,k)*(-0.5*p->DZN[KP]+0.5*fabs(dsz))-alpha(i,j,k);
            recheck=0.5*(nx(i,j,k)*p->DXN[KP]+ny(i,j,k)*p->DYN[JP]+nz(i,j,k)*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_m(i,j,k)-=Vol;
                V_w_p(i,j,k-1)-=Vol;
            }
            else if(nz(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_m(i,j,k)-=Vol;
                V_w_p(i,j,k-1)-=Vol;
            }
            
        }   
    }
}

void VOF_PLIC::advectWater_MACHO2D
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    // 2D
    if(sweep==0)
    {
        if(a->u(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->u(i,j,k)*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_p(i,j,k)+=Vol;
            V_w_m(i+1,j,k)+=Vol;
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double Vol;
            Vol=fabs(a->u(i-1,j,k))*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_m(i,j,k)-=Vol;
            V_w_p(i-1,j,k)-=Vol;
        }
    }
    else if(sweep==2)
    {
        if(a->w(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->w(i,j,k)*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_p(i,j,k)+=Vol;
            V_w_m(i,j,k+1)+=Vol;
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double Vol;
            Vol=fabs(a->w(i,j,k-1))*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_m(i,j,k)-=Vol;
            V_w_p(i,j,k-1)-=Vol;
        }
    }
}

void VOF_PLIC::advectPlane_Wang
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    if(sweep==0)
    {
        double QL,QR,UR,UL,nx_mod,r0_mod,nsum;
        UL=a->u(i-1,j,k);
        UR=a->u(i,j,k);
        QL=(UR*p->dt-UL*p->dt)*UL*p->dt/(2*p->DXN[IP])+UL*p->dt;
        QR=(UR*p->dt-UL*p->dt)*UR*p->dt/(2*p->DXN[IP])+UR*p->dt;
        nx_mod=nx(i,j,k)/(1.0+(QR-QL)/p->DXN[IP]);
        r0_mod=alpha(i,j,k)+nx_mod*QL;
        nsum=sqrt(nx_mod*nx_mod+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
        nx(i,j,k)=nx_mod/nsum;
        ny(i,j,k)=ny(i,j,k)/nsum;
        nz(i,j,k)=nz(i,j,k)/nsum;
        alpha(i,j,k)=r0_mod/nsum;
    }
    else if(sweep==1)
    {
        double QL,QR,UR,UL,ny_mod,r0_mod,nsum;
        UL=a->v(i,j-1,k);
        UR=a->v(i,j,k);
        QL=(UR*p->dt-UL*p->dt)*UL*p->dt/(2*p->DYN[JP])+UL*p->dt;
        QR=(UR*p->dt-UL*p->dt)*UR*p->dt/(2*p->DYN[JP])+UR*p->dt;
        ny_mod=ny(i,j,k)/(1.0+(QR-QL)/p->DYN[JP]);
        r0_mod=alpha(i,j,k)+ny_mod*QL;
        nsum=sqrt(ny_mod*ny_mod+nx(i,j,k)*nx(i,j,k)+nz(i,j,k)*nz(i,j,k));
        ny(i,j,k)=ny_mod/nsum;
        nx(i,j,k)=nx(i,j,k)/nsum;
        nz(i,j,k)=nz(i,j,k)/nsum;
        alpha(i,j,k)=r0_mod/nsum;
    }
    else
    {
        double QL,QR,UR,UL,nz_mod,r0_mod,nsum;
        UL=a->w(i,j,k-1);
        UR=a->w(i,j,k);
        QL=(UR*p->dt-UL*p->dt)*UL*p->dt/(2*p->DZN[KP])+UL*p->dt;
        QR=(UR*p->dt-UL*p->dt)*UR*p->dt/(2*p->DZN[KP])+UR*p->dt;
        nz_mod=nz(i,j,k)/(1.0+(QR-QL)/p->DZN[KP]);
        r0_mod=alpha(i,j,k)+nz_mod*QL;
        nsum=sqrt(nz_mod*nz_mod+nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k));
        nz(i,j,k)=nz_mod/nsum;
        nx(i,j,k)=nx(i,j,k)/nsum;
        ny(i,j,k)=ny(i,j,k)/nsum;
        alpha(i,j,k)=r0_mod/nsum;
    }
}

void VOF_PLIC::advectPlane_forBonnScheme
(
    fdm* a,
    lexer* p,
    int sweep
)
{

    if(sweep==0)
    {
        if(a->u(i,j,k)>0.0)
        {
            double dsx, scaledVol, Vol, r0x, recheck;
            dsx=a->u(i,j,k)*p->dt;
            r0x=nx(i,j,k)*(0.5*p->DXN[IP]-0.5*dsx)-alpha(i,j,k);
            recheck=0.5*(fabs(nx(i,j,k))*fabs(dsx)+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
                V_w_m(i+1,j,k)+=Vol;
            }
            else if(nx(i,j,k)<0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
                V_w_m(i+1,j,k)+=Vol;
            }
            
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double dsx, scaledVol, Vol, r0x,recheck;
            dsx=a->u(i-1,j,k)*p->dt;
            r0x=nx(i,j,k)*(-0.5*p->DXN[IP]+0.5*fabs(dsx))-alpha(i,j,k);
            recheck=0.5*(fabs(nx(i,j,k))*fabs(dsx)+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
                V_w_p(i-1,j,k)-=Vol;
            }
            else if(nx(i,j,k)>0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
                V_w_p(i-1,j,k)-=Vol;
            }
        }
    }
    else if(sweep==1)
    {
        if(a->v(i,j,k)>0.0)
        {
            double dsy,scaledVol, Vol, r0y,recheck;
            dsy=a->v(i,j,k)*p->dt;
            r0y=ny(i,j,k)*(0.5*p->DYN[JP]-0.5*dsy)-alpha(i,j,k);
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[IP]+fabs(ny(i,j,k))*fabs(dsy)+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0y);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],fabs(dsy),p->DZN[KP],r0y);
                Vol=scaledVol*p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
                V_w_m(i,j+1,k)+=Vol;
            }
            else if(ny(i,j,k)<0.0)
            {
                Vol=p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
                V_w_m(i,j+1,k)+=Vol;
            }
        }
        
        if(a->v(i,j-1,k)<0.0)
        {
            double dsy,scaledVol, Vol, r0y,recheck;
            dsy=a->v(i,j-1,k);
            r0y=ny(i,j,k)*(-0.5*p->DYN[JP]+0.5*fabs(dsy))-alpha(i,j,k);
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[IP]+fabs(ny(i,j,k))*fabs(dsy)+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0y);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],fabs(dsy),p->DZN[KP],r0y);
                Vol=scaledVol*p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
                V_w_p(i,j-1,k)-=Vol;
            }
            else if(ny(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
                V_w_p(i,j-1,k)-=Vol;
            }
        }
    }
    else
    {
        if(a->w(i,j,k)>0.0)
        {
            double dsz,scaledVol, Vol, r0z, recheck;
            dsz=a->w(i,j,k)*p->dt;
            r0z=nz(i,j,k)*(0.5*p->DZN[KP]-0.5*dsz)-alpha(i,j,k);
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[KP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_p(i,j,k)+=Vol;
                V_w_m(i,j,k+1)+=Vol;
            }
            else if(nz(i,j,k)<0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_p(i,j,k)+=Vol;
                V_w_m(i,j,k+1)+=Vol;
            }
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double dsz,scaledVol, Vol, r0z,recheck;
            dsz=a->w(i,j,k-1)*p->dt;
            r0z=nz(i,j,k)*(-0.5*p->DZN[KP]+0.5*fabs(dsz))-alpha(i,j,k);
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[KP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_m(i,j,k)-=Vol;
                V_w_p(i,j,k-1)-=Vol;
            }
            else if(nz(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_m(i,j,k)-=Vol;
                V_w_p(i,j,k-1)-=Vol;
            }
            
        }   
    }
}

void VOF_PLIC::advectWater_forBonnScheme
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    if(sweep==0)
    {
        if(a->u(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->u(i,j,k)*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_p(i,j,k)+=Vol;
            V_w_m(i+1,j,k)+=Vol;
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double Vol;
            Vol=fabs(a->u(i-1,j,k))*p->dt*p->DYN[JP]*p->DZN[KP];
            V_w_m(i,j,k)-=Vol;
            V_w_p(i-1,j,k)-=Vol;
        }
        
    }
    else if(sweep==1)
    {
        if(a->v(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->v(i,j,k)*p->dt*p->DXN[IP]*p->DZN[KP];
            V_w_p(i,j,k)+=Vol;
            V_w_m(i,j+1,k)+=Vol;
        }
        
        if(a->v(i,j-1,k)<0.0)
        {
            double Vol;
            Vol=fabs(a->v(i,j-1,k))*p->dt*p->DXN[IP]*p->DZN[KP];
            V_w_m(i,j,k)-=Vol;
            V_w_p(i,j-1,k)-=Vol;
        }
    }
    else
    {
        if(a->w(i,j,k)>0.0)
        {
            double Vol;
            Vol=a->w(i,j,k)*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_p(i,j,k)+=Vol;
            V_w_m(i,j,k+1)+=Vol;
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double Vol;
            Vol=fabs(a->w(i,j,k-1))*p->dt*p->DXN[IP]*p->DYN[JP];
            V_w_m(i,j,k)-=Vol;
            V_w_p(i,j,k-1)-=Vol;
        }
    }
}
    
    
        
    
    