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

void VOF_PLIC::advectPlane_forCOSMIC2D_RK
(
    fdm* a,
    lexer* p,
    int sweep,
    int inputdim,
    field& uvel,
    field& vvel,
    field& wvel
)
{
  
    if(sweep==0)
    {
        if(uvel(i,j,k)>0.0)
        {
            double dsx, scaledVol, Vol,r0x,recheck;
            if(p->F90==5)
                dsx=twoStepVel(p,a,uvel(i,j,k),uvel(i-1,j,k),p->DXN[IP]);
            else
                dsx=uvel(i,j,k)*p->dt;
            r0x=-(nx(i,j,k)*(0.5*p->DXN[IP]-0.5*dsx)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*fabs(dsx)+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0x);
            if(recheck>1E-20)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),dsx,p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*dsx*p->DYN[JP]*p->DZN[KP];
                if(Vol!=Vol)
                cout<<"plane u volnan in "<<i<<","<<j<<","<<k<<endl;
                switch(inputdim)
                {
                    case -1:
                                Vn_p(i,j,k)=Vol;
                                break;
                    case 0:
                                Vx_p(i,j,k)=Vol;
                                break;
                    case 2:
                                Vz_p(i,j,k)=Vol;
                                break;
                }
                
            }
            else if(nx(i,j,k)<0.0)
            {
                Vol=dsx*p->DYN[JP]*p->DZN[KP];
                switch(inputdim)
                {
                    case -1:
                                Vn_p(i,j,k)=Vol;
                                break;
                    case 0:
                                Vx_p(i,j,k)=Vol;
                                break;
                    case 2:
                                Vz_p(i,j,k)=Vol;
                                break;
                }
                
            }
        }
        
        if(uvel(i-1,j,k)<0.0)
        {
            double dsx, scaledVol, Vol, r0x, recheck;
            if(p->F90==5)
                dsx=fabs(twoStepVel(p,a,uvel(i-1,j,k),uvel(i,j,k),p->DXN[IP]));
            else
                dsx=fabs(uvel(i-1,j,k)*p->dt);
            r0x=-(nx(i,j,k)*(-0.5*p->DXN[IP]+0.5*dsx)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*fabs(dsx)+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*p->DZN[KP])-fabs(r0x);
            if(recheck>1E-20)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),dsx,p->DYN[JP],p->DZN[KP],r0x);
                Vol=scaledVol*dsx*p->DYN[JP]*p->DZN[KP];
                if(Vol!=Vol)
                    cout<<"plane u volnan in "<<i<<","<<j<<","<<k<<endl;
                switch(inputdim)
                {
                    case -1:
                                Vn_m(i,j,k)=-Vol;
                                break;
                    case 0:
                                Vx_m(i,j,k)=-Vol;
                                break;
                    case 2:
                                Vz_m(i,j,k)=-Vol;
                                break;
                }
            }
            else if(nx(i,j,k)>0.0)
            {
                Vol=dsx*p->DYN[JP]*p->DZN[KP];
                switch(inputdim)
                {
                    case -1:
                                Vn_m(i,j,k)=-Vol;
                                    break;
                    case 0:
                                Vx_m(i,j,k)=-Vol;
                                break;
                    case 2:
                                Vz_m(i,j,k)=-Vol;
                                break;
                }
                
            }
        }
        
    }
    else if(sweep==2)
    {
        if(wvel(i,j,k)>0.0)
        {
            double dsz, scaledVol, Vol, r0z, recheck;
            if(p->F90==5)
                dsz=twoStepVel(p,a,wvel(i,j,k),wvel(i,j,k-1),p->DZN[KP]);
            else
                dsz=wvel(i,j,k)*p->dt;
            r0z=-(nz(i,j,k)*(0.5*p->DZN[KP]-0.5*dsz)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*dsz)-fabs(r0z);
            if(recheck>1E-20)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],dsz,r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*dsz;
                if(Vol!=Vol)
                    cout<<"plane w volnan in "<<i<<","<<j<<","<<k<<endl;
                switch(inputdim)
                {
                    case -1:
                                Vn_p(i,j,k)=Vol;
                                break;
                    case 0:
                                Vx_p(i,j,k)=Vol;
                                break;
                    case 2:
                                Vz_p(i,j,k)=Vol;
                                break;
                }
            }
            else if(nz(i,j,k)<0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*dsz;
                switch(inputdim)
                {
                    case -1:
                                Vn_p(i,j,k)=Vol;
                                break;
                    case 0:
                                Vx_p(i,j,k)=Vol;
                                break;
                    case 2:
                                Vz_p(i,j,k)=Vol;
                                break;
                }
            }
        }
        
        if(wvel(i,j,k-1)<0.0)
        {
            double dsz, scaledVol, Vol, r0z, recheck;
            if(p->F90==5)
                dsz=fabs(twoStepVel(p,a,wvel(i,j,k-1),wvel(i,j,k),p->DZN[KP]));
            else
                dsz=fabs(wvel(i,j,k-1)*p->dt);
            r0z=-(nz(i,j,k)*(-0.5*p->DZN[KP]+0.5*dsz)-alpha(i,j,k));
            recheck=0.5*(fabs(nx(i,j,k))*p->DXN[IP]+fabs(ny(i,j,k))*p->DYN[JP]+fabs(nz(i,j,k))*dsz)-fabs(r0z);
            if(recheck>1E-20)
            {
                scaledVol=calculateVolume(nx(i,j,k),ny(i,j,k),nz(i,j,k),p->DXN[IP],p->DYN[JP],dsz,r0z);
                Vol=scaledVol*p->DXN[IP]*p->DYN[JP]*dsz;
                if(Vol!=Vol)
                    cout<<"plane w volnan in "<<i<<","<<j<<","<<k<<endl;
                switch(inputdim)
                {
                    case -1:
                                Vn_m(i,j,k)=-Vol;
                                break;
                    case 0:
                                Vx_m(i,j,k)=-Vol;
                                break;
                    case 2:
                                Vz_m(i,j,k)=-Vol;
                                break;
                }
            }
            else if(nz(i,j,k)>0.0)
            {
                Vol=p->DXN[IP]*p->DYN[JP]*dsz;
                switch(inputdim)
                {
                    case -1:
                                Vn_m(i,j,k)=-Vol;
                                break;
                    case 0:
                                Vx_m(i,j,k)=-Vol;
                                break;
                    case 2:
                                Vz_m(i,j,k)=-Vol;
                                break;
                }
            }
        }
        
    }
}

void VOF_PLIC::advectWater_forCOSMIC2D_RK
(
    fdm* a,
    lexer* p,
    int sweep,
    int inputdim,
    field& uvel,
    field& vvel,
    field& wvel
)
{
    if(sweep==0)
    {
        if(uvel(i,j,k)>0.0)
        {
            double dsx, Vol;
            if(p->F90==5)
                dsx=twoStepVel(p,a,uvel(i,j,k),uvel(i-1,j,k),p->DXN[IP]);
            else
                dsx=uvel(i,j,k)*p->dt;
            Vol=dsx*p->DYN[JP]*p->DZN[KP];
            if(Vol!=Vol)
                cout<<"water u volnan in "<<i<<","<<j<<","<<k<<endl;
            switch(inputdim)
            {
                case -1:
                            Vn_p(i,j,k)=Vol;
                            break;
                case 0:
                            Vx_p(i,j,k)=Vol;
                            break;
                case 2:
                            Vz_p(i,j,k)=Vol;
                            break;
            }
        }
        
        if(uvel(i-1,j,k)<0.0)
        {
            double dsx, Vol;
            if(p->F90==5)
                dsx=fabs(twoStepVel(p,a,uvel(i-1,j,k),uvel(i,j,k),p->DXN[IP]));
            else
                dsx=fabs(uvel(i-1,j,k)*p->dt);
            Vol=dsx*p->DYN[JP]*p->DZN[KP];
            if(Vol!=Vol)
                cout<<"water u volnan in "<<i<<","<<j<<","<<k<<endl;
            switch(inputdim)
            {
                case -1:
                            Vn_m(i,j,k)=-Vol;
                            break;
                case 0:
                            Vx_m(i,j,k)=-Vol;
                            break;
                case 2:
                            Vz_m(i,j,k)=-Vol;
                            break;
            }
        }
    }
    else if(sweep==2)
    {
        if(wvel(i,j,k)>0.0)
        {
            double dsz,Vol;
            if(p->F90==5)
                dsz=twoStepVel(p,a,wvel(i,j,k),wvel(i,j,k-1),p->DZN[KP]);
            else
                dsz=wvel(i,j,k)*p->dt;
            Vol=p->DXN[IP]*p->DYN[JP]*dsz;
            if(Vol!=Vol)
                cout<<"water w volnan in "<<i<<","<<j<<","<<k<<endl;
            switch(inputdim)
            {
                case -1:
                            Vn_p(i,j,k)=Vol;
                            break;
                case 0:
                            Vx_p(i,j,k)=Vol;
                            break;
                case 2:
                            Vz_p(i,j,k)=Vol;
                            break;
            }
        }
        
        if(wvel(i,j,k-1)<0.0)
        {
            double dsz,Vol;
            if(p->F90==5)
                dsz=fabs(twoStepVel(p,a,wvel(i,j,k-1),wvel(i,j,k),p->DZN[KP]));
            else
                dsz=fabs(wvel(i,j,k-1)*p->dt);
            Vol=p->DXN[IP]*p->DYN[JP]*dsz;
            if(Vol!=Vol)
                cout<<"water w volnan in "<<i<<","<<j<<","<<k<<endl;
            switch(inputdim)
            {
                case -1:
                            Vn_m(i,j,k)=-Vol;
                            break;
                case 0:
                            Vx_m(i,j,k)=-Vol;
                            break;
                case 2:
                            Vz_m(i,j,k)=-Vol;
                            break;
            }
        }
    }
}

double VOF_PLIC::twoStepVel(lexer* p, fdm* a,double uin1,double uin2,double deltax)
{
    double dx1, u12, dx2, dx_ret;
    dx1=0.5*p->dt*uin1;
    u12=(uin2-uin1)/deltax * dx1+uin1;
    dx2=0.5*p->dt*u12;
    dx_ret=dx1+dx2;
    
    return dx_ret;
}