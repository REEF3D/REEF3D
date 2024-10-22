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
    
void VOF_PLIC::advectPlane_NewWang
(
    fdm* a,
    lexer* p,
    int sweep
)
{
    if(sweep==0)
    {
        double xar,xal,xdr,xdl,xmedr,xmedl;
        double xdr_c,xdl_c,dlt_x;
        double u_xmedr,u_xmedl;
        int cellcnt_r,cellcnt_l;
        double Volsum=0.0;
        
        xar=0.5*p->DXN[IP];
        xal=-0.5*p->DXN[IP];
        xmedr=xar-0.5*p->dt*a->u(i,j,k);
        xmedl=xal-0.5*p->dt*a->u(i-1,j,k);
        u_xmedr=a->u(i,j,k)+(a->u(i+1,j,k)-a->u(i-1,j,k))/(p->DXN[IP]+p->DXN[IP1])*(-0.5*p->dt*a->u(i,j,k));
        u_xmedl=a->u(i-1,j,k)+(a->u(i,j,k)-a->u(i-2,j,k))/(p->DXN[IP]+p->DXN[IM1])*(-0.5*p->dt*a->u(i-1,j,k));
        xdr=xar-p->dt*u_xmedr;
        xdl=xal-p->dt*u_xmedl;
        
        if(xdr<=xar)
        {
            if(xdr<-0.5*p->DXN[IP])
            {
                if(xdr<(-0.5*p->DXN[IP]-p->DXN[IM1]))
                {
                    cellcnt_r=-2;
                    xdr_c=(p->DXP[IM1]+p->DXP[IM2])-fabs(xdr);
                }
                else
                {
                    cellcnt_r=-1;
                    xdr_c=p->DXP[IM1]-fabs(xdr);
                }
            }
            else
            {
                cellcnt_r=0;
                xdr_c=xdr;
            }
        }
        else
        {
            if(xdr>0.5*p->DXN[IP])
            {
                if(xdr>(0.5*p->DXN[IP]+p->DXN[IP1]))
                {
                    if(xdr>(0.5*p->DXN[IP]+p->DXN[IP1]+p->DXN[IP2]))
                    {
                        cellcnt_r=3;
                        xdr_c=-((p->DXP[IP]+p->DXP[IP1]+p->DXP[IP2])-fabs(xdr));
                    }
                    else
                    {
                        cellcnt_r=2;
                        xdr_c=-((p->DXP[IP]+p->DXP[IP1])-fabs(xdr));
                    }
                }
                else
                {
                    cellcnt_r=1;
                    xdr_c=-(p->DXP[IP]-fabs(xdr));
                }
            }
            else
            {
                cellcnt_r=0;
                xdr_c=xdr;
            }
        }
        
        if(xdl>=xal)
        {
            if(xdl>0.5*p->DXN[IP])
            {
                if(xdl>(0.5*p->DXN[IP]+p->DXN[IP1]))
                {
                    cellcnt_l=2;
                    xdl_c=-((p->DXP[IP]+p->DXP[IP1])-fabs(xdl));
                }
                else
                {
                    cellcnt_l=1;
                    xdl_c=-(p->DXP[IP]-fabs(xdl));
                }
            }
            else
            {
                cellcnt_l=0;
                xdl_c=xdl;
            }
        }
        else
        {
            if(xdl<(-0.5*p->DXN[IP]))
            {
                if(xdl<(-0.5*p->DXN[IP]-p->DXN[IM1]))
                {
                    if(xdl<(-0.5*p->DXN[IP]-p->DXN[IM1]-p->DXN[IM2]))
                    {
                        cellcnt_l=-3;
                        xdl_c=(p->DXP[IM1]-p->DXP[IM2]-p->DXP[IM3])-fabs(xdl);
                    }
                    else
                    {
                        cellcnt_l=-2;
                        xdl_c=(p->DXP[IM1]-p->DXP[IM2])-fabs(xdl);
                    }
                }
                else
                {
                    cellcnt_l=-1;
                    xdl_c=p->DXP[IM1]-fabs(xdl);
                }
            }
            else
            {
                cellcnt_l=0;
                xdl_c=xdl;
            }
        }
        
        for(int cellrun=cellcnt_l; cellrun<cellcnt_r; cellrun++)
        {   
           
            if(cellcnt_l==cellcnt_r)
                dlt_x=xdr_c-xdl_c;
            else if(cellrun==cellcnt_l)
            {
                switch(cellcnt_l)
                {
                    case -3:
                        dlt_x=0.5*p->DXN[IM3]-xdl_c;
                        break;
                    case -2:
                        dlt_x=0.5*p->DXN[IM2]-xdl_c;
                        break;
                    case -1:
                        dlt_x=0.5*p->DXN[IM1]-xdl_c;
                        break;
                    case 0:
                        dlt_x=0.5*p->DXN[IP]-xdl_c;
                        break;
                    case 1:
                        dlt_x=0.5*p->DXN[IP1]-xdl_c;
                        break;
                    case 2:
                        dlt_x=0.5*p->DXN[IP2]-xdl_c;
                        break;
                }
            }
            else if(cellrun==cellcnt_r)
            {
                switch(cellcnt_r)
                {
                    case -2:
                        dlt_x=0.5*p->DXN[IM2]+xdr_c;
                        break;
                    case -1:
                        dlt_x=0.5*p->DXN[IM1]+xdr_c;
                        break;
                    case 0:
                        dlt_x=0.5*p->DXN[IP]+xdr_c;
                        break;
                    case 1:
                        dlt_x=0.5*p->DXN[IP1]+xdr_c;
                        break;
                    case 2:
                        dlt_x=0.5*p->DXN[IP2]+xdr_c;
                        break;
                    case 3:
                        dlt_x=0.5*p->DXN[IP3]+xdr_c;
                        break;
                }
            }
            else
            {
                switch(cellrun)
                {   
                    case -3:
                        dlt_x=p->DXN[IM3];
                        break;
                    case -2:
                        dlt_x=p->DXN[IM2];
                        break;
                    case -1:
                        dlt_x=p->DXN[IM1];
                        break;
                    case 0:
                        dlt_x=p->DXN[IP];
                        break;
                    case 1:
                        dlt_x=p->DXN[IP1];
                        break;
                    case 2:
                        dlt_x=p->DXN[IP2];
                        break;
                    case 3:
                        dlt_x=p->DXN[IP3];
                        break;
                }
            }
            
            if(vofstep(i+cellrun,j,k)>0.9999)
                Volsum+=dlt_x*p->DYN[JP]*p->DZN[KP];
            else if(vofstep(i+cellrun,j,k)<0.0001)
                Volsum+=0.0;
            else
            {
                double r0xn,recheck,xn,cellVol,cellVol_s;
                if(cellcnt_l==cellcnt_r)
                    xn=0.5*(xdr_c+xdl_c);
                else if(cellrun==cellcnt_l)
                    xn=xdl_c+0.5*dlt_x;
                else if(cellrun==cellcnt_r)
                    xn=xdr_c-0.5*dlt_x;
                else
                    xn=0.0;
                    
                r0xn=-(nx(i+cellrun,j,k)*xn-alpha(i+cellrun,j,k));
                recheck=0.5*(fabs(nx(i+cellrun,j,k))*dlt_x+fabs(ny(i+cellrun,j,k))*p->DYN[JP]+nz(i+cellrun,j,k)*p->DZN[KP])-fabs(r0xn);
                if(recheck>0.0)
                {
                    cellVol_s=calculateVolume(nx(i+cellrun,j,k),ny(i+cellrun,j,k),nz(i+cellrun,j,k),dlt_x,p->DYN[JP],p->DZN[KP],r0xn);
                    cellVol=cellVol_s*dlt_x*p->DYN[JP]*p->DZN[KP];
                    Volsum+=cellVol;
                }
            }
        }
        V_w_p(i,j,k)=Volsum/(fabs(xdr-xdl)*p->DYN[JP]*p->DZN[KP]);
    }
    
    else if(sweep==1)
    {
        double yar,yal,ydr,ydl,ymedr,ymedl;
        double ydr_c,ydl_c,dlt_y;
        double v_ymedr,v_ymedl;
        int cellcnt_r,cellcnt_l;
        double Volsum=0.0;
        
        yar=0.5*p->DYN[JP];
        yal=-0.5*p->DYN[JP];
        ymedr=yar-0.5*p->dt*a->v(i,j,k);
        ymedl=yal-0.5*p->dt*a->v(i,j-1,k);
        v_ymedr=a->v(i,j,k)+(a->v(i,j+1,k)-a->v(i,j-1,k))/(p->DYN[JP]+p->DYN[JP1])*(-0.5*p->dt*a->v(i,j,k));
        v_ymedl=a->v(i,j-1,k)+(a->v(i,j,k)-a->v(i,j-2,k))/(p->DYN[JP]+p->DYN[JM1])*(-0.5*p->dt*a->v(i,j-1,k));
        ydr=yar-p->dt*v_ymedr;
        ydl=yal-p->dt*v_ymedl;
        
        if(ydr<=yar)
        {
            if(ydr<-0.5*p->DYN[JP])
            {
                if(ydr<(-0.5*p->DYN[JP]-p->DYN[JM1]))
                {
                    cellcnt_r=-2;
                    ydr_c=(p->DYP[JM1]+p->DYP[JM2])-fabs(ydr);
                }
                else
                {
                    cellcnt_r=-1;
                    ydr_c=p->DYP[JM1]-fabs(ydr);
                }
            }
            else
            {
                cellcnt_r=0;
                ydr_c=ydr;
            }
        }
        else
        {
            if(ydr>0.5*p->DYN[JP])
            {
                if(ydr>(0.5*p->DYN[JP]+p->DYN[JP1]))
                {
                    if(ydr>(0.5*p->DYN[JP]+p->DYN[JP1]+p->DYN[JP2]))
                    {
                        cellcnt_r=3;
                        ydr_c=-((p->DYP[JP]+p->DYP[JP1]+p->DYP[JP2])-fabs(ydr));
                    }
                    else
                    {
                        cellcnt_r=2;
                        ydr_c=-((p->DYP[JP]+p->DYP[JP1])-fabs(ydr));
                    }
                }
                else
                {
                    cellcnt_r=1;
                    ydr_c=-(p->DYP[JP]-fabs(ydr));
                }
            }
            else
            {
                cellcnt_r=0;
                ydr_c=ydr;
            }
        }
        
        if(ydl>=yal)
        {
            if(ydl>0.5*p->DYN[JP])
            {
                if(ydl>(0.5*p->DYN[JP]+p->DYN[JP1]))
                {
                    cellcnt_l=2;
                    ydl_c=-((p->DYP[JP]+p->DYP[JP1])-fabs(ydl));
                }
                else
                {
                    cellcnt_l=1;
                    ydl_c=-(p->DYP[JP]-fabs(ydl));
                }
            }
            else
            {
                cellcnt_l=0;
                ydl_c=ydl;
            }
        }
        else
        {
            if(ydl<(-0.5*p->DYN[JP]))
            {
                if(ydl<(-0.5*p->DYN[JP]-p->DYN[JM1]))
                {
                    if(ydl<(-0.5*p->DYN[JP]-p->DYN[JM1]-p->DYN[JM2]))
                    {
                        cellcnt_l=-3;
                        ydl_c=(p->DYP[JM1]-p->DYP[JM2]-p->DYP[JM3])-fabs(ydl);
                    }
                    else
                    {
                        cellcnt_l=-2;
                        ydl_c=(p->DYP[JM1]-p->DYP[JM2])-fabs(ydl);
                    }
                }
                else
                {
                    cellcnt_l=-1;
                    ydl_c=p->DYP[JM1]-fabs(ydl);
                }
            }
            else
            {
                cellcnt_l=0;
                ydl_c=ydl;
            }
        }
        
        for(int cellrun=cellcnt_l; cellrun<cellcnt_r; cellrun++)
        {   
           
            if(cellcnt_l==cellcnt_r)
                dlt_y=ydr_c-ydl_c;
            else if(cellrun==cellcnt_l)
            {
                switch(cellcnt_l)
                {
                    case -3:
                        dlt_y=0.5*p->DYN[JM3]-ydl_c;
                        break;
                    case -2:
                        dlt_y=0.5*p->DYN[JM2]-ydl_c;
                        break;
                    case -1:
                        dlt_y=0.5*p->DYN[JM1]-ydl_c;
                        break;
                    case 0:
                        dlt_y=0.5*p->DYN[JP]-ydl_c;
                        break;
                    case 1:
                        dlt_y=0.5*p->DYN[JP1]-ydl_c;
                        break;
                    case 2:
                        dlt_y=0.5*p->DYN[JP2]-ydl_c;
                        break;
                }
            }
            else if(cellrun==cellcnt_r)
            {
                switch(cellcnt_r)
                {
                    case -2:
                        dlt_y=0.5*p->DYN[JM2]+ydr_c;
                        break;
                    case -1:
                        dlt_y=0.5*p->DYN[JM1]+ydr_c;
                        break;
                    case 0:
                        dlt_y=0.5*p->DYN[JP]+ydr_c;
                        break;
                    case 1:
                        dlt_y=0.5*p->DYN[JP1]+ydr_c;
                        break;
                    case 2:
                        dlt_y=0.5*p->DYN[JP2]+ydr_c;
                        break;
                    case 3:
                        dlt_y=0.5*p->DYN[JP3]+ydr_c;
                        break;
                }
            }
            else
            {
                switch(cellrun)
                {   
                    case -3:
                        dlt_y=p->DYN[JM3];
                        break;
                    case -2:
                        dlt_y=p->DYN[JM2];
                        break;
                    case -1:
                        dlt_y=p->DYN[JM1];
                        break;
                    case 0:
                        dlt_y=p->DYN[JP];
                        break;
                    case 1:
                        dlt_y=p->DYN[JP1];
                        break;
                    case 2:
                        dlt_y=p->DYN[JP2];
                        break;
                    case 3:
                        dlt_y=p->DYN[JP3];
                        break;
                }
            }
            
            if(vofstep(i,j+cellrun,k)>0.9999)
                Volsum+=p->DXN[IP]*dlt_y*p->DZN[KP];
            else if(vofstep(i,j+cellrun,k)<0.0001)
                Volsum+=0.0;
            else
            {
                double r0yn,recheck,yn,cellVol,cellVol_s;
                if(cellcnt_l==cellcnt_r)
                    yn=0.5*(ydr_c+ydl_c);
                else if(cellrun==cellcnt_l)
                    yn=ydl_c+0.5*dlt_y;
                else if(cellrun==cellcnt_r)
                    yn=ydr_c-0.5*dlt_y;
                else
                    yn=0.0;
                    
                r0yn=-(ny(i,j+cellrun,k)*yn-alpha(i,j+cellrun,k));
                recheck=0.5*(fabs(nx(i,j+cellrun,k))*p->DXN[IP]+fabs(ny(i,j+cellrun,k))*dlt_y+nz(i,j+cellrun,k)*p->DZN[KP])-fabs(r0yn);
                if(recheck>0.0)
                {
                    cellVol_s=calculateVolume(nx(i,j+cellrun,k),ny(i,j+cellrun,k),nz(i,j+cellrun,k),p->DXN[IP],dlt_y,p->DZN[KP],r0yn);
                    cellVol=cellVol_s*p->DXN[IP]*dlt_y*p->DZN[KP];
                    Volsum+=cellVol;
                }
            }
        }
        V_w_p(i,j,k)=Volsum/(p->DXN[IP]*fabs(ydr-ydl)*p->DZN[KP]);
    }
    
    else
    {
        double zar,zal,zdr,zdl,zmedr,zmedl;
        double zdr_c,zdl_c,dlt_z;
        double w_zmedr,w_zmedl;
        int cellcnt_r,cellcnt_l;
        double Volsum=0.0;
        
        zar=0.5*p->DZN[KP];
        zal=-0.5*p->DZN[KP];
        zmedr=zar-0.5*p->dt*a->w(i,j,k);
        zmedl=zal-0.5*p->dt*a->w(i,j,k-1);
        w_zmedr=a->w(i,j,k)+(a->w(i,j,k+1)-a->w(i,j,k-1))/(p->DZN[KP]+p->DZN[KP1])*(-0.5*p->dt*a->w(i,j,k));
        w_zmedl=a->w(i,j,k-1)+(a->w(i,j,k)-a->w(i,j,k-2))/(p->DZN[KP]+p->DZN[KM1])*(-0.5*p->dt*a->w(i,j,k-1));
        zdr=zar-p->dt*w_zmedr;
        zdl=zal-p->dt*w_zmedl;
        
        if(zdr<=zar)
        {
            if(zdr<-0.5*p->DZN[KP])
            {
                if(zdr<(-0.5*p->DZN[KP]-p->DZN[KM1]))
                {
                    cellcnt_r=-2;
                    zdr_c=(p->DZP[KM1]+p->DZP[KM2])-fabs(zdr);
                }
                else
                {
                    cellcnt_r=-1;
                    zdr_c=p->DZP[KM1]-fabs(zdr);
                }
            }
            else
            {
                cellcnt_r=0;
                zdr_c=zdr;
            }
        }
        else
        {
            if(zdr>0.5*p->DZN[KP])
            {
                if(zdr>(0.5*p->DZN[KP]+p->DZN[KP1]))
                {
                    if(zdr>(0.5*p->DZN[KP]+p->DZN[KP1]+p->DZN[KP2]))
                    {
                        cellcnt_r=3;
                        zdr_c=-((p->DZP[KP]+p->DZP[KP1]+p->DZP[KP2])-fabs(zdr));
                    }
                    else
                    {
                        cellcnt_r=2;
                        zdr_c=-((p->DZP[KP]+p->DZP[KP1])-fabs(zdr));
                    }
                }
                else
                {
                    cellcnt_r=1;
                    zdr_c=-(p->DZP[KP]-fabs(zdr));
                }
            }
            else
            {
                cellcnt_r=0;
                zdr_c=zdr;
            }
        }
        
        if(zdl>=zal)
        {
            if(zdl>0.5*p->DZN[KP])
            {
                if(zdl>(0.5*p->DZN[KP]+p->DZN[KP1]))
                {
                    cellcnt_l=2;
                    zdl_c=-((p->DZP[KP]+p->DZP[KP1])-fabs(zdl));
                }
                else
                {
                    cellcnt_l=1;
                    zdl_c=-(p->DZP[KP]-fabs(zdl));
                }
            }
            else
            {
                cellcnt_l=0;
                zdl_c=zdl;
            }
        }
        else
        {
            if(zdl<(-0.5*p->DZN[KP]))
            {
                if(zdl<(-0.5*p->DZN[KP]-p->DZN[KM1]))
                {
                    if(zdl<(-0.5*p->DZN[KP]-p->DZN[KM1]-p->DZN[KM2]))
                    {
                        cellcnt_l=-3;
                        zdl_c=(p->DZP[KM1]-p->DZP[KM2]-p->DZP[KM3])-fabs(zdl);
                    }
                    else
                    {
                        cellcnt_l=-2;
                        zdl_c=(p->DZP[KM1]-p->DZP[KM2])-fabs(zdl);
                    }
                }
                else
                {
                    cellcnt_l=-1;
                    zdl_c=p->DZP[KM1]-fabs(zdl);
                }
            }
            else
            {
                cellcnt_l=0;
                zdl_c=zdl;
            }
        }
        
        for(int cellrun=cellcnt_l; cellrun<cellcnt_r; cellrun++)
        {   
           
            if(cellcnt_l==cellcnt_r)
                dlt_z=zdr_c-zdl_c;
            else if(cellrun==cellcnt_l)
            {
                switch(cellcnt_l)
                {
                    case -3:
                        dlt_z=0.5*p->DZN[KM3]-zdl_c;
                        break;
                    case -2:
                        dlt_z=0.5*p->DZN[KM2]-zdl_c;
                        break;
                    case -1:
                        dlt_z=0.5*p->DZN[KM1]-zdl_c;
                        break;
                    case 0:
                        dlt_z=0.5*p->DZN[KP]-zdl_c;
                        break;
                    case 1:
                        dlt_z=0.5*p->DZN[KP1]-zdl_c;
                        break;
                    case 2:
                        dlt_z=0.5*p->DZN[KP2]-zdl_c;
                        break;
                }
            }
            else if(cellrun==cellcnt_r)
            {
                switch(cellcnt_r)
                {
                    case -2:
                        dlt_z=0.5*p->DZN[KM2]+zdr_c;
                        break;
                    case -1:
                        dlt_z=0.5*p->DZN[KM1]+zdr_c;
                        break;
                    case 0:
                        dlt_z=0.5*p->DZN[KP]+zdr_c;
                        break;
                    case 1:
                        dlt_z=0.5*p->DZN[KP1]+zdr_c;
                        break;
                    case 2:
                        dlt_z=0.5*p->DZN[KP2]+zdr_c;
                        break;
                    case 3:
                        dlt_z=0.5*p->DZN[KP3]+zdr_c;
                        break;
                }
            }
            else
            {
                switch(cellrun)
                {   
                    case -3:
                        dlt_z=p->DZN[KM3];
                        break;
                    case -2:
                        dlt_z=p->DZN[KM2];
                        break;
                    case -1:
                        dlt_z=p->DZN[KM1];
                        break;
                    case 0:
                        dlt_z=p->DZN[KP];
                        break;
                    case 1:
                        dlt_z=p->DZN[KP1];
                        break;
                    case 2:
                        dlt_z=p->DZN[KP2];
                        break;
                    case 3:
                        dlt_z=p->DZN[KP3];
                        break;
                }
            }
            
            if(vofstep(i,j,k+cellrun)>0.9999)
                Volsum+=p->DXN[IP]*p->DYN[JP]*dlt_z;
            else if(vofstep(i,j+cellrun,k)<0.0001)
                Volsum+=0.0;
            else
            {
                double r0zn,recheck,zn,cellVol,cellVol_s;
                if(cellcnt_l==cellcnt_r)
                    zn=0.5*(zdr_c+zdl_c);
                else if(cellrun==cellcnt_l)
                    zn=zdl_c+0.5*dlt_z;
                else if(cellrun==cellcnt_r)
                    zn=zdr_c-0.5*dlt_z;
                else
                    zn=0.0;
                    
                r0zn=-(nz(i,j,k+cellrun)*zn-alpha(i,j,k+cellrun));
                recheck=0.5*(fabs(nx(i,j,k+cellrun))*p->DXN[IP]+fabs(ny(i,j,k+cellrun))*p->DYN[JP]+nz(i,j,k+cellrun)*dlt_z)-fabs(r0zn);
                if(recheck>0.0)
                {
                    cellVol_s=calculateVolume(nx(i,j,k+cellrun),ny(i,j,k+cellrun),nz(i,j,k+cellrun),p->DXN[IP],p->DYN[JP],dlt_z,r0zn);
                    cellVol=cellVol_s*p->DXN[IP]*p->DYN[JP]*dlt_z;
                    Volsum+=cellVol;
                }
            }
        }
        V_w_p(i,j,k)=Volsum/(p->DXN[IP]*p->DYN[JP]*fabs(zdr-zdl));
    }
}
        
    
    