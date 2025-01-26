/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"fnpf_fsf_update.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"slice.h"

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e20)

void fnpf_fsf_update::velcalc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *f)
{/*
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
    */
    FLOOP
    {
    // U
    if(k<p->knoz)
    c->U[FIJK] = (-c->Fi[FIp2JK] + 8.0*c->Fi[FIp1JK] - 8.0*c->Fi[FIm1JK] + c->Fi[FIm2JK])/(-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2])
    
                + p->sigx[FIJK]*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));
                
                
    if(k==p->knoz)
    c->U[FIJK] = (c->Fi[FIp1JK]-c->Fi[FIm1JK])/(p->DXP[IP]+p->DXP[IM1])
    
                + p->sigx[FIJK]*((c->Fi[FIJK]-c->Fi[FIJKm1])/(p->DZN[KP]));
     
    // V           
    if(k<p->knoz)
    c->V[FIJK] = (-c->Fi[FIJp2K] + 8.0*c->Fi[FIJp1K] - 8.0*c->Fi[FIJm1K] + c->Fi[FIJm2K])/(-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2])
                
                + p->sigy[FIJK]*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));
                
    if(k==p->knoz)
    c->V[FIJK] = (c->Fi[FIJp1K]-c->Fi[FIJm1K])/(p->DYP[JP]+p->DYP[JM1])
                
                + p->sigy[FIJK]*((c->Fi[FIJK]-c->Fi[FIJKm1])/(p->DZN[KP]));
    
    // W
    c->W[FIJK] = ((-c->Fi[FIJKp2] + 8.0*c->Fi[FIJKp1] - 8.0*c->Fi[FIJKm1] + c->Fi[FIJKm2])/(-p->ZN[KP2] + 8.0*p->ZN[KP1] - 8.0*p->ZN[KM1] + p->ZN[KM2]))*p->sigz[IJ];
    }

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

