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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_weno3.h"
#include"lexer.h"
#include"vec.h"
#include"field.h"
#include"fnpf_discrete_weights.h"

fnpf_weno3::fnpf_weno3(lexer* p) :  ddweno3_f_nug(p)
{
    p->Darray(ckz,p->knoz+1+4*marge,5);
    
    fnpf_discrete_weights dw(p);

    dw.ck_weights(p, ckz, p->ZN, p->knoz+1, 1, 4, 6);
}

fnpf_weno3::~fnpf_weno3()
{
}

double fnpf_weno3::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    grad=0.0;
    
    if(0.5*(ivel1+ivel2)>0.0)
    grad=ddwenox(f,1.0);
    
    if(0.5*(ivel1+ivel2)<0.0)
    grad=ddwenox(f,-1.0);
    
    return grad;
}

double fnpf_weno3::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    grad=0.0;
    
    if(0.5*(jvel1+jvel2)>0.0)
    grad=ddwenoy(f,1.0);
    
    if(0.5*(jvel1+jvel2)<0.0)
    grad=ddwenoy(f,-1.0);
    
    return grad;
}

double fnpf_weno3::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;
    

    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3] && p->flag4[IJKm4]>0 && p->flag4[IJKm5])
    {
        if(i+p->origin_i>0)
        grad = (-(49.0/20.0)*f(i,j,k+1) + 6.0*f(i,j,k) - 7.5*f(i,j,k-1) + (20.0/3.0)*f(i,j,k-2) - (15.0/4.0)*f(i,j,k-3) + (6.0/5.0)*f(i,j,k-4) - (1.0/6.0)*f(i,j,k-5))
          /(-(49.0/20.0)*p->ZP[KP1] + 6.0*p->ZP[KP] - 7.5*p->ZP[KM1] + (20.0/3.0)*p->ZP[KM2] - (15.0/4.0)*p->ZP[KM3] + (6.0/5.0)*p->ZP[KM4] - (1.0/6.0)*p->ZP[KM5]);
              
        if(i+p->origin_i==0)
        grad = (-(49.0/20.0)*f(i,j,k) + 6.0*f(i,j,k-1) - 7.5*f(i,j,k-2) + (20.0/3.0)*f(i,j,k-3) - (15.0/4.0)*f(i,j,k-4) + (6.0/5.0)*f(i,j,k-5) - (1.0/6.0)*f(i,j,k-6))
          /(-(49.0/20.0)*p->ZP[KP] + 6.0*p->ZP[KM1] - 7.5*p->ZP[KM2] + (20.0/3.0)*p->ZP[KM3] - (15.0/4.0)*p->ZP[KM4] + (6.0/5.0)*p->ZP[KM5] - (1.0/6.0)*p->ZP[KM6]);
              
        //cout<<" return 6"<<endl;
            
        return grad;
    }

    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3]>0)
    {
        if(i+p->origin_i>0)
        grad = (-(25.0/12.0)*f(i,j,k+1) + 4.0*f(i,j,k) - 3.0*f(i,j,k-1) + (4.0/3.0)*f(i,j,k-2) - 0.25*f(i,j,k-3))
              /(-(25.0/12.0)*p->ZP[KP1] + 4.0*p->ZP[KP] - 3.0*p->ZP[KM1] + (4.0/3.0)*p->ZP[KM2] - 0.25*p->ZP[KM3]);
              
        if(i+p->origin_i==0)
        grad = (-(25.0/12.0)*f(i,j,k) + 4.0*f(i,j,k-1) - 3.0*f(i,j,k-2) + (4.0/3.0)*f(i,j,k-3) - 0.25*f(i,j,k-4))
              /(-(25.0/12.0)*p->ZP[KP] + 4.0*p->ZP[KM1] - 3.0*p->ZP[KM2] + (4.0/3.0)*p->ZP[KM3] - 0.25*p->ZP[KM4]);
              
        //cout<<" return 4"<<endl;
            
        return grad;
    }
    
    else
    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0)
    {
        if(i+p->origin_i>0)
        grad = (-1.5*f(i,j,k+1) + 2.0*f(i,j,k) - 0.5*f(i,j,k-1))/(-1.5*p->ZP[KP1] + 2.0*p->ZP[KP] - 0.5*p->ZP[KM1]);
              
        if(i+p->origin_i==0)
        grad = (-1.5*f(i,j,k) + 2.0*f(i,j,k-1) - 0.5*f(i,j,k-2))/(-1.5*p->ZP[KP] + 2.0*p->ZP[KM1] - 0.5*p->ZP[KM2]);
            
        //cout<<" return 2"<<endl;    
        
        return grad;
    }
    
    else
    {
        if(i+p->origin_i>0)
        grad = (f(i,j,k+1) - f(i,j,k))/(p->ZP[KP]);
              
        if(i+p->origin_i==0)
        grad = (f(i,j,k) - f(i,j,k-1))/(p->ZP[KM1]);
            
        //cout<<" return 1"<<endl;    
            
        return grad;
    }
}

double fnpf_weno3::sx(lexer *p, slice &f, double ivel)
{
    grad=0.0;
    
    if(ivel>0.0)
    grad=dswenox(f,1.0);
    
    if(ivel<0.0)
    grad=dswenox(f,-1.0);
    
    return grad;
}

double fnpf_weno3::sy(lexer *p, slice &f, double jvel)
{
    grad=0.0;
    
    if(jvel>0.0)
    grad=dswenoy(f,1.0);
    
    if(jvel<0.0)
    grad=dswenoy(f,-1.0);
    
    return grad;   
}

double fnpf_weno3::sz(lexer *p, double *f)
{
    grad = (ckz[p->knoz+marge][4]*f[FIJK] + ckz[p->knoz+marge][3]*f[FIJKm1] + ckz[p->knoz+marge][2]*f[FIJKm2] 
          + ckz[p->knoz+marge][1]*f[FIJKm3] + ckz[p->knoz+marge][0]*f[FIJKm4]);
          
    return grad;   
}
