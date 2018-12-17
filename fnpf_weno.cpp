/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_weno.h"
#include"lexer.h"
#include"vec.h"

fnpf_weno::fnpf_weno(lexer* p) :  ddweno_f_nug(p)
{
}

fnpf_weno::~fnpf_weno()
{
}

double fnpf_weno::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    grad=0.0;
    
    if(0.5*(ivel1+ivel2)>0.0)
    grad=ddwenox(f,1.0);
    
    if(0.5*(ivel1+ivel2)<0.0)
    grad=ddwenox(f,-1.0);
    
    return grad;
}

double fnpf_weno::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    grad=0.0;
    
    if(0.5*(jvel1+jvel2)>0.0)
    grad=ddwenoy(f,1.0);
    
    if(0.5*(jvel1+jvel2)<0.0)
    grad=ddwenoy(f,-1.0);
    
    return grad;
}

double fnpf_weno::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;
    
    if(0.5*(kvel1+kvel2)>0.0)
    grad=ddwenoz(f,1.0);
    
    if(0.5*(kvel1+kvel2)<0.0)
    grad=ddwenoz(f,-1.0);
    
    return grad;
}

double fnpf_weno::sx(lexer *p, slice &f, double ivel)
{
    grad=0.0;
    
    if(ivel>0.0)
    grad=dswenox(f,1.0);
    
    if(ivel<0.0)
    grad=dswenox(f,-1.0);
    
    return grad;
}

double fnpf_weno::sy(lexer *p, slice &f, double jvel)
{
    grad=0.0;
    
    if(jvel>0.0)
    grad=dswenoy(f,1.0);
    
    if(jvel<0.0)
    grad=dswenoy(f,-1.0);
    
    return grad;   
}

double fnpf_weno::sz(lexer *p, double *f)
{
    return (-(25.0/12.0)*f[FIJK] + 4.0*f[FIJKm1] - 3.0*f[FIJKm2] + (4.0/3.0)*f[FIJKm3] - 0.25*f[FIJKm4])
          /(-(25.0/12.0)*p->ZN[KP] + 4.0*p->ZN[KM1] - 3.0*p->ZN[KM2] + (4.0/3.0)*p->ZN[KM3] - 0.25*p->ZN[KM4]); 
    /*
    return (-(25.0/12.0)*f[FIJKp1] + 4.0*f[FIJK] - 3.0*f[FIJKm1] + (4.0/3.0)*f[FIJKm2] - 0.25*f[FIJKm3])
          /(-(25.0/12.0)*p->ZN[KP1] + 4.0*p->ZN[KP] - 3.0*p->ZN[KM1] + (4.0/3.0)*p->ZN[KM2] - 0.25*p->ZN[KM3]);
           */ 
}
