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


void VOF_PLIC::advectPlane_forBonnScheme
(
    fdm* a,
    lexer* p,
    int sweep
)
{

    if(sweep==0)
    {
        double usx; /*
        if(nx(i,j,k)>0.0)
            usx=a->u(i-1,j,k);
        else
            usx=a->u(i,j,k);*/
        if(a->u(i,j,k)>0.0)
        {
            double  dsx, scaledVol, Vol, r0x, recheck;
            usx=a->u(i,j,k);
            //usx=a->u(i,j,k)-0.5*a->u(i,j,k)*p->dt*(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
            dsx=usx*p->dt;
            r0x=-(nx(i,j,k)*(0.5*p->DXN[IP]-0.5*dsx)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*fabs(dsx)+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
            }
            else if(nx(i,j,k)<0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
            }
            
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double  dsx, scaledVol, Vol, r0x,recheck;
            usx=a->u(i-1,j,k);
            //usx=a->u(i-1,j,k)-0.5*a->u(i-1,j,k)*p->dt*(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
            dsx=usx*p->dt;
            r0x=-(nx(i,j,k)*(-0.5*p->DXN[IP]+0.5*fabs(dsx))-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*fabs(dsx)+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0x);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),fabs(dsx),p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
            }
            else if(nx(i,j,k)>0.0)
            {
                Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
            }
        }
    }
    else if(sweep==1)
    {
        if(a->v(i,j,k)>0.0)
        {
            double vsy, dsy,scaledVol, Vol, r0y,recheck;
//            vsy=a->v(i,j,k)-0.5*a->v(i,j,k)*p->dt*(a->v(i,j,k)-a->v(i,j-1,k))/p->DYN[JP];
            dsy=a->v(i,j,k)*p->dt;
            r0y=-(ny(i,j,k)*(0.5*p->DYN[JP]-0.5*dsy)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[IP]+fabs(ny(i,j,k))*fabs(dsy)+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0y);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],fabs(dsy),p->DZN[KP],r0y);
                Vol=scaledVol*p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
            }
            else if(ny(i,j,k)<0.0)
            {
                Vol=p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_p(i,j,k)+=Vol;
            }
        }
        
        if(a->v(i,j-1,k)<0.0)
        {
            double vsy, dsy,scaledVol, Vol, r0y,recheck;
          //  vsy=a->v(i,j-1,k)-0.5*a->v(i,j-1,k)*p->dt*(a->v(i,j,k)-a->v(i,j-1,k))/p->DYN[JP];
            dsy=a->v(i,j-1,k)*p->dt;
            r0y=-(ny(i,j,k)*(-0.5*p->DYN[JP]+0.5*fabs(dsy))-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[IP]+fabs(ny(i,j,k))*fabs(dsy)+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0y);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],fabs(dsy),p->DZN[KP],r0y);
                Vol=scaledVol*p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
            }
            else if(ny(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*fabs(dsy)*p->DZN[KP];
                V_w_m(i,j,k)-=Vol;
            }
        }
    }
    else
    {
        double wsz; /*
        if(nz(i,j,k)>=0.0)
            wsz=a->w(i,j,k-1);
        else
            wsz=a->w(i,j,k);*/
            
        if(a->w(i,j,k)>0.0)
        {
            double dsz,scaledVol, Vol, r0z, recheck;
            wsz=a->w(i,j,k);
            //wsz=a->w(i,j,k)-0.5*a->w(i,j,k)*p->dt*(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
            dsz=wsz*p->dt;
            r0z=-(nz(i,j,k)*(0.5*p->DZN[KP]-0.5*dsz)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[KP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_p(i,j,k)+=Vol;
            }
            else if(nz(i,j,k)<0.0)
            {
               Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_p(i,j,k)+=Vol;
            }
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double  dsz,scaledVol, Vol, r0z,recheck;
            wsz=a->w(i,j,k-1);
            //wsz=a->w(i,j,k-1)-0.5*a->w(i,j,k-1)*p->dt*(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
            dsz=wsz*p->dt;
            r0z=-(nz(i,j,k)*(-0.5*p->DZN[KP]+0.5*fabs(dsz))-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[KP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*fabs(dsz))-fabs(r0z);
            if(recheck>0.0)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],fabs(dsz),r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_m(i,j,k)-=Vol;

           }
           else if(nz(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*fabs(dsz);
                V_w_m(i,j,k)-=Vol;
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
            double usx, dsx, Vol;
           // usx=a->u(i,j,k)-0.5*a->u(i,j,k)*p->dt*(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
            dsx=a->u(i,j,k)*p->dt;
            Vol=dsx*p->DYN[JP]*p->DZN[KP];
            V_w_p(i,j,k)+=Vol;
        }
        if(a->u(i-1,j,k)<0.0)
        {
            double usx, dsx, Vol;
           // usx=a->u(i-1,j,k)-0.5*a->u(i-1,j,k)*p->dt*(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
            dsx=a->u(i-1,j,k)*p->dt;
            Vol=fabs(dsx)*p->DYN[JP]*p->DZN[KP];
            V_w_m(i,j,k)-=Vol;
        }
        
    }
    else if(sweep==1)
    {
        if(a->v(i,j,k)>0.0)
        {
            double vsy, dsy, Vol;
            //vsy=a->v(i,j,k)-0.5*a->v(i,j,k)*p->dt*(a->v(i,j,k)-a->v(i,j-1,k))/p->DYN[JP];
            dsy=a->v(i,j,k)*p->dt;
            Vol=dsy*p->DXN[IP]*p->DZN[KP];
            V_w_p(i,j,k)+=Vol;
        }
        
        if(a->v(i,j-1,k)<0.0)
        {
            double vsy, dsy, Vol;
           // vsy=a->v(i,j-1,k)-0.5*a->v(i,j-1,k)*p->dt*(a->v(i,j,k)-a->v(i,j-1,k))/p->DYN[JP];
            dsy=a->v(i,j-1,k)*p->dt;
            Vol=fabs(dsy)*p->DXN[IP]*p->DZN[KP];
            V_w_m(i,j,k)-=Vol;
        }
    }
    else
    {
        if(a->w(i,j,k)>0.0)
        {
            double wsz,dsz, Vol;
           // wsz=a->w(i,j,k)-0.5*a->w(i,j,k)*p->dt*(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
            dsz=a->w(i,j,k)*p->dt;
            Vol=dsz*p->DXN[IP]*p->DYN[JP];
            V_w_p(i,j,k)+=Vol;
        }
        if(a->w(i,j,k-1)<0.0)
        {
            double wsz, dsz, Vol;
           // wsz=a->w(i,j,k-1)-0.5*a->w(i,j,k-1)*p->dt*(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
            dsz=a->w(i,j,k-1)*p->dt;
            Vol=fabs(dsz)*p->DXN[IP]*p->DYN[JP];
            V_w_m(i,j,k)-=Vol;
        }
    }
}
    
    
        
    
    