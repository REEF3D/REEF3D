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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_ioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    
    // Uio
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->Uio_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    k=patch[qq]->gcb[n][2];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        u(i-1,j,k) =  patch[qq]->Uio;
        u(i-2,j,k) =  patch[qq]->Uio;
        u(i-3,j,k) =  patch[qq]->Uio;
        
        v(i-1,j,k) =  0.0;
        v(i-2,j,k) =  0.0;
        v(i-3,j,k) =  0.0;
        
        w(i-1,j,k) =  0.0;
        w(i-2,j,k) =  0.0;
        w(i-3,j,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        u(i,j+1,k) =  0.0;
        u(i,j+2,k) =  0.0;
        u(i,j+3,k) =  0.0;
        
        v(i,j,k)   =  patch[qq]->Uio;
        v(i,j+1,k) =  patch[qq]->Uio;
        v(i,j+2,k) =  patch[qq]->Uio;
        
        w(i,j+1,k) =  0.0;
        w(i,j+2,k) =  0.0;
        w(i,j+3,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        u(i,j-1,k) =  0.0;
        u(i,j-2,k) =  0.0;
        u(i,j-3,k) =  0.0;
        
        v(i,j-1,k) =  patch[qq]->Uio;
        v(i,j-2,k) =  patch[qq]->Uio;
        v(i,j-3,k) =  patch[qq]->Uio;
        
        w(i,j-1,k) =  0.0;
        w(i,j-2,k) =  0.0;
        w(i,j-3,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        u(i,j,k)   =  patch[qq]->Uio;
        u(i+1,j,k) =  patch[qq]->Uio;
        u(i+2,j,k) =  patch[qq]->Uio;
        
        v(i+1,j,k) =  0.0;
        v(i+2,j,k) =  0.0;
        v(i+3,j,k) =  0.0;
        
        w(i+1,j,k) =  0.0;
        w(i+2,j,k) =  0.0;
        w(i+3,j,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==5)
        {
        u(i,j,k-1) =  0.0;
        u(i,j,k-2) =  0.0;
        u(i,j,k-3) =  0.0;
        
        v(i,j,k-1) =  0.0;
        v(i,j,k-2) =  0.0;
        v(i,j,k-3) =  0.0;
        
        w(i,j,k-1) =  patch[qq]->Uio;
        w(i,j,k-2) =  patch[qq]->Uio;
        w(i,j,k-3) =  patch[qq]->Uio;
        }
        
        if(patch[qq]->gcb[n][3]==6)
        {
        u(i,j,k+1) =  0.0;
        u(i,j,k+2) =  0.0;
        u(i,j,k+3) =  0.0;
        
        v(i,j,k+1) =  0.0;
        v(i,j,k+2) =  0.0;
        v(i,j,k+3) =  0.0;
        
        w(i,j,k)   =  patch[qq]->Uio;
        w(i,j,k+1) =  patch[qq]->Uio;
        w(i,j,k+2) =  patch[qq]->Uio;
        }
    
    }
    
    // Discharge
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->Q_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    k=patch[qq]->gcb[n][2];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        u(i-1,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        u(i-2,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        u(i-3,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        v(i-1,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        v(i-2,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        v(i-3,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        w(i-1,j,k) =  0.0;
        w(i-2,j,k) =  0.0;
        w(i-3,j,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        u(i,j+1,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        u(i,j+2,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        u(i,j+3,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        v(i,j,k)   =  patch[qq]->Uq*patch[qq]->sinalpha;
        v(i,j+1,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        v(i,j+2,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        w(i,j+1,k) =  0.0;
        w(i,j+2,k) =  0.0;
        w(i,j+3,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        u(i,j-1,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        u(i,j-2,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        u(i,j-3,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        v(i,j-1,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        v(i,j-2,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        v(i,j-3,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        w(i,j-1,k) =  0.0;
        w(i,j-2,k) =  0.0;
        w(i,j-3,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        u(i,j,k)   =  patch[qq]->Uq*patch[qq]->sinalpha;
        u(i+1,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        u(i+2,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        v(i+1,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        v(i+2,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        v(i+3,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        w(i+1,j,k) =  0.0;
        w(i+2,j,k) =  0.0;
        w(i+3,j,k) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==5)
        {
        u(i,j,k-1) =  0.0;
        u(i,j,k-2) =  0.0;
        u(i,j,k-3) =  0.0;
        
        v(i,j,k-1) =  0.0;
        v(i,j,k-2) =  0.0;
        v(i,j,k-3) =  0.0;
        
        w(i,j,k-1) =  patch[qq]->Uq;
        w(i,j,k-2) =  patch[qq]->Uq;
        w(i,j,k-3) =  patch[qq]->Uq;
        }
        
        if(patch[qq]->gcb[n][3]==6)
        {
        u(i,j,k+1) =  0.0;
        u(i,j,k+2) =  0.0;
        u(i,j,k+3) =  0.0;
        
        v(i,j,k+1) =  0.0;
        v(i,j,k+2) =  0.0;
        v(i,j,k+3) =  0.0;
        
        w(i,j,k)   =  patch[qq]->Uq;
        w(i,j,k+1) =  patch[qq]->Uq;
        w(i,j,k+2) =  patch[qq]->Uq;
        }
    
    }
    
    
    
     // Velocity components
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->velcomp_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    k=patch[qq]->gcb[n][2];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        u(i-1,j,k) =  patch[qq]->U;
        u(i-2,j,k) =  patch[qq]->U;
        u(i-3,j,k) =  patch[qq]->U;
        
        v(i,j,k)   =  patch[qq]->V;
        v(i-1,j,k) =  patch[qq]->V;
        v(i-2,j,k) =  patch[qq]->V;
        v(i-3,j,k) =  patch[qq]->V;
        
        w(i,j,k)   =  patch[qq]->W;
        w(i-1,j,k) =  patch[qq]->W;
        w(i-2,j,k) =  patch[qq]->W;
        w(i-3,j,k) =  patch[qq]->W;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        u(i,j,k)   =  patch[qq]->U;
        u(i,j+1,k) =  patch[qq]->U;
        u(i,j+2,k) =  patch[qq]->U;
        u(i,j+3,k) =  patch[qq]->U;
        
        v(i,j,k)   =  patch[qq]->V;
        v(i,j+1,k) =  patch[qq]->V;
        v(i,j+2,k) =  patch[qq]->V;
        
        w(i,j,k)   =  patch[qq]->W;
        w(i,j+1,k) =  patch[qq]->W;
        w(i,j+2,k) =  patch[qq]->W;
        w(i,j+3,k) =  patch[qq]->W;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        u(i,j,k)   =  patch[qq]->U;
        u(i,j-1,k) =  patch[qq]->U;
        u(i,j-2,k) =  patch[qq]->U;
        u(i,j-3,k) =  patch[qq]->U;
        
        v(i,j-1,k) =  patch[qq]->V;
        v(i,j-2,k) =  patch[qq]->V;
        v(i,j-3,k) =  patch[qq]->V;
        
        w(i,j,k)   =  patch[qq]->W;
        w(i,j-1,k) =  patch[qq]->W;
        w(i,j-2,k) =  patch[qq]->W;
        w(i,j-3,k) =  patch[qq]->W;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        u(i,j,k)   =  patch[qq]->U;
        u(i+1,j,k) =  patch[qq]->U;
        u(i+2,j,k) =  patch[qq]->U;
        
        v(i,j,k)   =  patch[qq]->V;
        v(i+1,j,k) =  patch[qq]->V;
        v(i+2,j,k) =  patch[qq]->V;
        v(i+3,j,k) =  patch[qq]->V;
        
        w(i,j,k)   =  patch[qq]->W;
        w(i+1,j,k) =  patch[qq]->W;
        w(i+2,j,k) =  patch[qq]->W;
        w(i+3,j,k) =  patch[qq]->W;
        }
        
        if(patch[qq]->gcb[n][3]==5)
        {
        u(i,j,k)   =  patch[qq]->U;
        u(i,j,k-1) =  patch[qq]->U;
        u(i,j,k-2) =  patch[qq]->U;
        u(i,j,k-3) =  patch[qq]->U;
        
        v(i,j,k)   =  patch[qq]->V;
        v(i,j,k-1) =  patch[qq]->V;
        v(i,j,k-2) =  patch[qq]->V;
        v(i,j,k-3) =  patch[qq]->V;
        
        w(i,j,k-1) =  patch[qq]->W;
        w(i,j,k-2) =  patch[qq]->W;
        w(i,j,k-3) =  patch[qq]->W;
        }
        
        if(patch[qq]->gcb[n][3]==6)
        {
        u(i,j,k)   =  patch[qq]->U;
        u(i,j,k+1) =  patch[qq]->U;
        u(i,j,k+2) =  patch[qq]->U;
        u(i,j,k+3) =  patch[qq]->U;
        
        v(i,j,k)   =  patch[qq]->V;
        v(i,j,k+1) =  patch[qq]->V;
        v(i,j,k+2) =  patch[qq]->V;
        v(i,j,k+3) =  patch[qq]->V;
        
        w(i,j,k)   =  patch[qq]->W;
        w(i,j,k+1) =  patch[qq]->W;
        w(i,j,k+2) =  patch[qq]->W;
        }
    
    }


} 

void patchBC::patchBC_rkioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    // Velocity components
    for(qq=0;qq<obj_count;++qq)
    //if(patch[qq]->velcomp_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    k=patch[qq]->gcb[n][2];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        u(i-1,j,k) =  a->u(i-1,j,k);
        u(i-2,j,k) =  a->u(i-2,j,k);
        u(i-3,j,k) =  a->u(i-3,j,k);
        
        v(i-1,j,k) =  a->v(i-1,j,k);
        v(i-2,j,k) =  a->v(i-2,j,k);
        v(i-3,j,k) =  a->v(i-3,j,k);
        
        w(i-1,j,k) =  a->w(i-1,j,k);
        w(i-2,j,k) =  a->w(i-2,j,k);
        w(i-3,j,k) =  a->w(i-3,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        u(i,j+1,k) =  a->u(i,j+1,k);
        u(i,j+2,k) =  a->u(i,j+2,k);
        u(i,j+3,k) =  a->u(i,j+3,k);
        
        v(i,j,k)   =  a->v(i,j,k);
        v(i,j+1,k) =  a->v(i,j+1,k);
        v(i,j+2,k) =  a->v(i,j+2,k);
        
        w(i,j+1,k) =  a->w(i,j+1,k);
        w(i,j+2,k) =  a->w(i,j+2,k);
        w(i,j+3,k) =  a->w(i,j+3,k);
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        u(i,j-1,k) =  a->u(i,j-1,k);
        u(i,j-2,k) =  a->u(i,j-2,k);
        u(i,j-3,k) =  a->u(i,j-3,k);
        
        v(i,j-1,k) =  a->v(i,j-1,k);
        v(i,j-2,k) =  a->v(i,j-2,k);
        v(i,j-3,k) =  a->v(i,j-3,k);
        
        w(i,j-1,k) =  a->w(i,j-1,k);
        w(i,j-2,k) =  a->w(i,j-2,k);
        w(i,j-3,k) =  a->w(i,j-3,k);
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        u(i,j,k)   =  a->u(i,j,k);
        u(i+1,j,k) =  a->u(i+1,j,k);
        u(i+2,j,k) =  a->u(i+2,j,k);
        
        v(i+1,j,k) =  a->v(i+1,j,k);
        v(i+2,j,k) =  a->v(i+2,j,k);
        v(i+3,j,k) =  a->v(i+3,j,k);
        
        w(i+1,j,k) =  a->w(i+1,j,k);
        w(i+2,j,k) =  a->w(i+2,j,k);
        w(i+3,j,k) =  a->w(i+3,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==5)
        {
        u(i,j,k-1) =  a->u(i,j,k-1);
        u(i,j,k-2) =  a->u(i,j,k-2);
        u(i,j,k-3) =  a->u(i,j,k-3);
        
        v(i,j,k-1) =  a->v(i,j,k-1);
        v(i,j,k-2) =  a->v(i,j,k-2);
        v(i,j,k-3) =  a->v(i,j,k-3);
        
        w(i,j,k-1) =  a->w(i,j,k-1);
        w(i,j,k-2) =  a->w(i,j,k-2);
        w(i,j,k-3) =  a->w(i,j,k-3);
        }
        
        if(patch[qq]->gcb[n][3]==6)
        {
        u(i,j,k+1) =  a->u(i,j,k+1);
        u(i,j,k+2) =  a->u(i,j,k+2);
        u(i,j,k+3) =  a->u(i,j,k+3);
        
        v(i,j,k+1) =  a->v(i,j,k+1);
        v(i,j,k+2) =  a->v(i,j,k+2);
        v(i,j,k+3) =  a->v(i,j,k+3);
        
        w(i,j,k)   =  a->w(i,j,k);
        w(i,j,k+1) =  a->w(i,j,k+1);
        w(i,j,k+2) =  a->w(i,j,k+2);
        }
    
    }

}

void patchBC::patchBC_ioflow2D(lexer *p, ghostcell*, slice&, slice&, slice&, slice&)
{
    cout<<p->mpirank<<" patchBC_2D: inflow 3D "<<endl;
}

void patchBC::patchBC_rkioflow2D(lexer *p, ghostcell*, slice&, slice&, slice&, slice&)
{
}


