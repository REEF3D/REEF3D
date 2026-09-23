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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include "fnpf_weno5.h"
#include "lexer.h"
#include "vec.h"
#include "field.h"
#include "fnpf_discrete_weights.h"

fnpf_weno5::fnpf_weno5(lexer *p) : ddweno_f_nug(p)
{
    p->Darray(ckz,p->knoz+1+4*marge,5);

    fnpf_discrete_weights dw(p);

    dw.ck_weights(p, ckz, p->ZN, p->knoz+1, 1, 4, 6);
}

fnpf_weno5::~fnpf_weno5()
{
}

double fnpf_weno5::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    if(0.5*(ivel1+ivel2)>0.0)
        return ddwenox(f,1.0);
    else if(0.5*(ivel1+ivel2)<0.0)
        return ddwenox(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    if(0.5*(jvel1+jvel2)>0.0)
        return ddwenoy(f,1.0);
    else if(0.5*(jvel1+jvel2)<0.0)
        return ddwenoy(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;

    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3] && p->flag4[IJKm4]>0 && p->flag4[IJKm5])
    {
        if(i+p->origin_i>0)
            grad = (-(49.0/20.0)*f(i,j,k+1) + 6.0*f(i,j,k) - 7.5*f(i,j,k-1) + (20.0/3.0)*f(i,j,k-2) - (15.0/4.0)*f(i,j,k-3) + (6.0/5.0)*f(i,j,k-4) - (1.0/6.0)*f(i,j,k-5))
            /(-(49.0/20.0)*p->ZP[KP1] + 6.0*p->ZP[KP] - 7.5*p->ZP[KM1] + (20.0/3.0)*p->ZP[KM2] - (15.0/4.0)*p->ZP[KM3] + (6.0/5.0)*p->ZP[KM4] - (1.0/6.0)*p->ZP[KM5]);
        else if(i+p->origin_i==0)
            grad = (-(49.0/20.0)*f(i,j,k) + 6.0*f(i,j,k-1) - 7.5*f(i,j,k-2) + (20.0/3.0)*f(i,j,k-3) - (15.0/4.0)*f(i,j,k-4) + (6.0/5.0)*f(i,j,k-5) - (1.0/6.0)*f(i,j,k-6))
            /(-(49.0/20.0)*p->ZP[KP] + 6.0*p->ZP[KM1] - 7.5*p->ZP[KM2] + (20.0/3.0)*p->ZP[KM3] - (15.0/4.0)*p->ZP[KM4] + (6.0/5.0)*p->ZP[KM5] - (1.0/6.0)*p->ZP[KM6]);

        return grad;
    }
    else if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3]>0)
    {
        if(i+p->origin_i>0)
            grad = (-(25.0/12.0)*f(i,j,k+1) + 4.0*f(i,j,k) - 3.0*f(i,j,k-1) + (4.0/3.0)*f(i,j,k-2) - 0.25*f(i,j,k-3))
                /(-(25.0/12.0)*p->ZP[KP1] + 4.0*p->ZP[KP] - 3.0*p->ZP[KM1] + (4.0/3.0)*p->ZP[KM2] - 0.25*p->ZP[KM3]);
        else if(i+p->origin_i==0)
            grad = (-(25.0/12.0)*f(i,j,k) + 4.0*f(i,j,k-1) - 3.0*f(i,j,k-2) + (4.0/3.0)*f(i,j,k-3) - 0.25*f(i,j,k-4))
                /(-(25.0/12.0)*p->ZP[KP] + 4.0*p->ZP[KM1] - 3.0*p->ZP[KM2] + (4.0/3.0)*p->ZP[KM3] - 0.25*p->ZP[KM4]);

        return grad;
    }
    else if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0)
    {
        if(i+p->origin_i>0)
            grad = (-1.5*f(i,j,k+1) + 2.0*f(i,j,k) - 0.5*f(i,j,k-1))/(-1.5*p->ZP[KP1] + 2.0*p->ZP[KP] - 0.5*p->ZP[KM1]);
        else if(i+p->origin_i==0)
            grad = (-1.5*f(i,j,k) + 2.0*f(i,j,k-1) - 0.5*f(i,j,k-2))/(-1.5*p->ZP[KP] + 2.0*p->ZP[KM1] - 0.5*p->ZP[KM2]);

        return grad;
    }
    else
    {
        if(i+p->origin_i>0)
            grad = (f(i,j,k+1) - f(i,j,k))/(p->ZP[KP]);
        else if(i+p->origin_i==0)
            grad = (f(i,j,k) - f(i,j,k-1))/(p->ZP[KM1]);

        return grad;
    }
}

double fnpf_weno5::sx(lexer *p, slice &f, double ivel)
{
    if(ivel>0.0)
        return dswenox(f,1.0);
    else if(ivel<0.0)
        return dswenox(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5::sy(lexer *p, slice &f, double jvel)
{
    if(jvel>0.0)
        return dswenoy(f,1.0);
    else if(jvel<0.0)
        return dswenoy(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5::sz(lexer *p, double *f)
{
    if(p->flag7[FIJK]>0 && p->flag7[FIJKm1]>0 && p->flag7[FIJKm2]>0 && p->flag7[FIJKm3]>0)
    {
        return (-(25.0/12.0)*f[FIJK] + 4.0*f[FIJKm1] - 3.0*f[FIJKm2] + (4.0/3.0)*f[FIJKm3] - 0.25*f[FIJKm4])
              /(-(25.0/12.0)*p->ZN[KP] + 4.0*p->ZN[KM1] - 3.0*p->ZN[KM2] + (4.0/3.0)*p->ZN[KM3] - 0.25*p->ZN[KM4]);
    }
    else if(p->flag7[FIJK]>0 && p->flag7[FIJKm1]>0)
    {
        return (-1.5*f[FIJK] + 2.0*f[FIJKm1] - 0.5*f[FIJKm2])/(-1.5*p->ZN[KP] + 2.0*p->ZN[KM1] - 0.5*p->ZN[KM2]);
    }
    else
    {
        return (f[FIJK] - f[FIJKm1])/(p->ZN[KM1]);
    }
}
