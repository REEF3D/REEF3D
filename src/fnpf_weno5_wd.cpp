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

#include"fnpf_weno5_wd.h"
#include"lexer.h"
#include"vec.h"
#include"fnpf_discrete_weights.h"

fnpf_weno5_wd::fnpf_weno5_wd(lexer *p,fdm_fnpf *c) :  fnpf_ddweno_f_nug(p,c)
{
    p->Darray(ckz,p->knoz+1+4*marge,5);
    
    fnpf_discrete_weights dw(p);

    dw.ck_weights(p, ckz, p->ZN, p->knoz+1, 1, 4, 6);
}

fnpf_weno5_wd::~fnpf_weno5_wd()
{
}

double fnpf_weno5_wd::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    grad=0.0;
    
    if(0.5*(ivel1+ivel2)>0.0)
    grad=ddwenox(f,1.0);
    
    if(0.5*(ivel1+ivel2)<0.0)
    grad=ddwenox(f,-1.0);
    
    return grad;
}

double fnpf_weno5_wd::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    grad=0.0;
    
    if(0.5*(jvel1+jvel2)>0.0)
    grad=ddwenoy(f,1.0);
    
    if(0.5*(jvel1+jvel2)<0.0)
    grad=ddwenoy(f,-1.0);
    
    return grad;
}

double fnpf_weno5_wd::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;
    
    if(0.5*(kvel1+kvel2)>0.0)
    grad=ddwenoz(f,1.0);
    
    if(0.5*(kvel1+kvel2)<0.0)
    grad=ddwenoz(f,-1.0);
    
    return grad;
}

double fnpf_weno5_wd::sx(lexer *p, slice &f, double ivel)
{
    grad=0.0;
    
    if(ivel>0.0)
    grad=dswenox(f,1.0);
    
    if(ivel<0.0)
    grad=dswenox(f,-1.0);
    
    return grad;
}

double fnpf_weno5_wd::sy(lexer *p, slice &f, double jvel)
{
    grad=0.0;
    
    if(jvel>0.0)
    grad=dswenoy(f,1.0);
    
    if(jvel<0.0)
    grad=dswenoy(f,-1.0);
    
    return grad;   
}

double fnpf_weno5_wd::sz(lexer *p, double *f)
{
    grad = (ckz[p->knoz+marge][4]*f[FIJK] + ckz[p->knoz+marge][3]*f[FIJKm1] + ckz[p->knoz+marge][2]*f[FIJKm2] 
        + ckz[p->knoz+marge][1]*f[FIJKm3] + ckz[p->knoz+marge][0]*f[FIJKm4]);
          
    //grad = (-(49.0/20.0)*f[FIJK] + 6.0*f[FIJKm1] - 7.5*f[FIJKm2] + (20.0/3.0)*f[FIJKm3] - 3.75*f[FIJKm4] + (6.0/5.0)*f[FIJKm5] - (1.0/6.0)*f[FIJKm6])
    //      /(-(49.0/20.0)*p->ZN[KP] + 6.0*p->ZN[KM1] - 7.5*p->ZN[KM2] + (20.0/3.0)*p->ZN[KM3] - 3.75*p->ZN[KM4] + (6.0/5.0)*p->ZN[KM5] - (1.0/6.0)*p->ZN[KM6]);
   
    
    return grad; 
    
}
