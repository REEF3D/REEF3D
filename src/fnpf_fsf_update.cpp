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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_fsf_update.h"#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"onephase.h"
#include"slice.h"

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e20)

fnpf_fsf_update::fnpf_fsf_update(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
    gcval_u = 10;
    gcval_v = 11;
    gcval_w = 12;
}

fnpf_fsf_update::~fnpf_fsf_update()
{
    
}

void fnpf_fsf_update::fsfupdate(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, onephase *poneph, slice &eta)
{
}

void fnpf_fsf_update::etaloc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void fnpf_fsf_update::etaloc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // find k location for eta
    SLICELOOP4
    c->etaloc(i,j) = p->knoz;
    
    pgc->gcsl_start4int(p,c->etaloc,50);
}

void fnpf_fsf_update::fsfbc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Fifsf, double *Fi)
{
    FFILOOP4
    {
        Fi[FIJK]   = Fifsf(i,j);
        Fi[FIJKp1] = Fifsf(i,j);
        Fi[FIJKp2] = Fifsf(i,j);  
        Fi[FIJKp3] = Fifsf(i,j);
    }
    
    /*
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        k=p->knoz;
        Fi[FIm1JKp1] = Fifsf(i,j);
    }*/
}

void fnpf_fsf_update::fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Fifsf, field &Fi)
{
}

void fnpf_fsf_update::fsfepol(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, field &Fi)
{
}

void fnpf_fsf_update::velcalc(lexer *p, fdm_fnpf *c, ghostcell *pgc, field &f)
{
}

void fnpf_fsf_update::velcalc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{
    FLOOP
    {
    // U
    if(k<p->knoz)
    c->U[FIJK] = (c->Fi[FIp1JK]-c->Fi[FIm1JK])/(p->DXP[IP]+p->DXP[IM1])
    
                + p->sigx[FIJK]*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));
                
    if(k==p->knoz)
    c->U[FIJK] = (c->Fi[FIp1JK]-c->Fi[FIm1JK])/(p->DXP[IP]+p->DXP[IM1])
    
                + p->sigx[FIJK]*((c->Fi[FIJK]-c->Fi[FIJKm1])/(p->DZN[KP]));
     
    // V           
    if(k<p->knoz)
    c->V[FIJK] = (c->Fi[FIJp1K]-c->Fi[FIJm1K])/(p->DYP[JP]+p->DYP[JM1])
                
                + p->sigy[FIJK]*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));
                
    if(k==p->knoz)
    c->V[FIJK] = (c->Fi[FIJp1K]-c->Fi[FIJm1K])/(p->DYP[JP]+p->DYP[JM1])
                
                + p->sigy[FIJK]*((c->Fi[FIJK]-c->Fi[FIJKm1])/(p->DZN[KP]));
    
    // W
    c->W[FIJK] = ((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZP[KP]+p->DZP[KM1]))*p->sigz[IJ];
    }
    
    //if(p->A343>=1)
    FLOOP
    {
        if(p->wet[Im1J]==0 || p->wet[Ip1J]==0 || p->wet[IJm1]==0 || p->wet[IJp1]==0 
        || p->wet[Im1Jm1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Ip1Jp1]==0)
        {
        
        c->U[FIJK]=0.0;
        c->V[FIJK]=0.0;
        c->W[FIJK]=0.0;
        }
    }

    /*
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        FKLOOP
        FPCHECK
        {
        c->U[FIJK] = c->Uin[FIm1JK];

        }
    }*/
        
    FFILOOP4
    c->W[FIJK] = c->Fz(i,j);
    
    
    
    int gcval=210;
    
    pgc->start7V(p,c->U,c->bc,gcval);
    pgc->start7V(p,c->V,c->bc,gcval);
    pgc->start7V(p,c->W,c->bc,gcval);
    
    
    
    // test: kfsfbc
    double val;
    
    
    SLICELOOP4
    {
    c->eta_n(i,j) = c->eta(i,j);
    }
    pgc->gcsl_start4(p,c->eta_n,1);   
    
}

