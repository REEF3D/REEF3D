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

#include"fnpf_sg_fsfbc.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"field4.h"
#include"discrete.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"onephase.h"
#include"fnpf_voiddisc.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"fnpf_cds6.h"
#include"fnpf_weno.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"

fnpf_sg_fsfbc::fnpf_sg_fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc) : Fx(p),Fy(p),Ex(p),Ey(p),
                                                                Exx(p),Eyy(p),Bx(p),By(p),Bxx(p),Byy(p)
{    
    if(p->A311==0)
    pdisc = new fnpf_voiddisc(p);
    
    if(p->A311==2)
    pdisc = new fnpf_cds2(p);
    
    if(p->A311==3)
    pdisc = new fnpf_cds4(p);
    
    if(p->A311==4)
    pdisc = new fnpf_weno(p);
    
    if(p->A311==6)
    pdisc = new fnpf_cds6(p);
    
    if(p->A312==2)
    pddx = new fnpf_ddx_cds2(p);
    
    if(p->A312==3)
    pddx = new fnpf_ddx_cds4(p);
}

fnpf_sg_fsfbc::~fnpf_sg_fsfbc()
{
}

void fnpf_sg_fsfbc::fsfdisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    c->WL(i,j) = eta(i,j) + p->wd - c->bed(i,j);
    
    // fi
    FFILOOP4
    {
    ivel1 = (Fifsf(i,j) - Fifsf(i-1,j))/p->DXP[IM1];
    ivel2 = (Fifsf(i+1,j) - Fifsf(i,j))/p->DXP[IP];
    
    jvel1 = (Fifsf(i,j) - Fifsf(i,j-1))/p->DYP[JM1];
    jvel2 = (Fifsf(i,j+1) - Fifsf(i,j))/p->DYP[JP];
    
    
    Fx(i,j) = pdisc->sx(p,Fifsf,ivel1,ivel2);
    Fy(i,j) = pdisc->sy(p,Fifsf,jvel1,jvel2);
    
    Ex(i,j) = pdisc->sx(p,eta,ivel1,ivel2);
    Ey(i,j) = pdisc->sy(p,eta,jvel1,jvel2);
    
    
    Exx(i,j) = pddx->sxx(p,eta);
    Eyy(i,j) = pddx->syy(p,eta);
    
    Bx(i,j) = pdisc->sx(p,c->depth,ivel1,ivel2);
    By(i,j) = pdisc->sy(p,c->depth,jvel1,jvel2);
    
    Bxx(i,j) = pddx->sxx(p,c->depth);
    Byy(i,j) = pddx->syy(p,c->depth);
    }
    
}

void fnpf_sg_fsfbc::fsfwvel(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    FFILOOP4
    {
    kvel =  p->sigz[FIJK]*(c->Fi[FIJK] - c->Fi[FIJKm1])/p->DZP[KM1];

    c->Fz(i,j) = p->sigz[FIJK]*pdisc->sz(p,c->Fi,kvel,kvel);
    }
}

void fnpf_sg_fsfbc::kfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    SLICELOOP4
    c->K(i,j) =  - Fx(i,j)*Ex(i,j) - Fy(i,j)*Ey(i,j) 
    
                 + c->Fz(i,j)*(1.0 + pow(Ex(i,j),2.0) + pow(Ey(i,j),2.0));
}

void fnpf_sg_fsfbc::dfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta)
{
    SLICELOOP4
    c->K(i,j) =  - 0.5*pow(Fx(i,j),2.0) - 0.5*pow(Fy(i,j),2.0) 
    
                 + 0.5*pow(c->Fz(i,j),2.0)*(1.0 + pow(Ex(i,j),2.0) + pow(Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j);
}


