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

#include "fnpf_weno5_wd.h"
#include "lexer.h"
#include "vec.h"
#include "fnpf_discrete_weights.h"

fnpf_weno5_wd::fnpf_weno5_wd(lexer *p) : fnpf_ddweno_f_nug(p), szw(p)
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
    if(0.5*(ivel1+ivel2)>0.0)
        return ddwenox(f,1.0);
    else if(0.5*(ivel1+ivel2)<0.0)
        return ddwenox(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5_wd::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    if(0.5*(jvel1+jvel2)>0.0)
        return ddwenoy(f,1.0);
    else if(0.5*(jvel1+jvel2)<0.0)
        return ddwenoy(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5_wd::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    if(0.5*(kvel1+kvel2)>0.0)
        return ddwenoz(f,1.0);
    else if(0.5*(kvel1+kvel2)<0.0)
        return ddwenoz(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5_wd::sx(lexer *p, slice &f, double ivel)
{
    if(ivel>0.0)
        return dswenox(f,1.0);
    else if(ivel<0.0)
        return dswenox(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5_wd::sy(lexer *p, slice &f, double jvel)
{
    if(jvel>0.0)
        return dswenoy(f,1.0);
    else if(jvel<0.0)
        return dswenoy(f,-1.0);
    else
        return 0.0;
}

double fnpf_weno5_wd::sz(lexer *p, double *f)
{
    return szw.eval(p,f,FIJK,5);
}
