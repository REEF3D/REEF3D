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

#include"patchBC_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC_2D::patchBC_ioflow2D(lexer *p, ghostcell *pgc, slice &P, slice &Q, slice &eta, slice &bed)
{
    // Uio
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->Uio_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        P(i,j) =  patch[qq]->Uio;
        P(i-1,j) =  patch[qq]->Uio;
        P(i-2,j) =  patch[qq]->Uio;
        P(i-3,j) =  patch[qq]->Uio;
        
        Q(i-1,j) =  0.0;
        Q(i-2,j) =  0.0;
        Q(i-3,j) =  0.0;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        P(i,j+1) =  0.0;
        P(i,j+2) =  0.0;
        P(i,j+3) =  0.0;
        
        Q(i,j)   =  patch[qq]->Uio;
        Q(i,j+1) =  patch[qq]->Uio;
        Q(i,j+2) =  patch[qq]->Uio;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        P(i,j-1) =  0.0;
        P(i,j-2) =  0.0;
        P(i,j-3) =  0.0;
        
        Q(i,j-1) =  patch[qq]->Uio;
        Q(i,j-2) =  patch[qq]->Uio;
        Q(i,j-3) =  patch[qq]->Uio;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        P(i,j)   =  patch[qq]->Uio;
        P(i+1,j) =  patch[qq]->Uio;
        P(i+2,j) =  patch[qq]->Uio;
        
        Q(i+1,j) =  0.0;
        Q(i+2,j) =  0.0;
        Q(i+3,j) =  0.0;
        }
    }
    
    // Uq discharge + Q hydrograph
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->Q_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        P(i-1,j) =  patch[qq]->Uq*patch[qq]->sinalpha;
        P(i-2,j) =  patch[qq]->Uq*patch[qq]->sinalpha;
        P(i-3,j) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        Q(i-1,j) =  patch[qq]->Uq*patch[qq]->cosalpha;
        Q(i-2,j) =  patch[qq]->Uq*patch[qq]->cosalpha;
        Q(i-3,j) =  patch[qq]->Uq*patch[qq]->cosalpha;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        P(i,j+1) =  patch[qq]->Uq*patch[qq]->cosalpha;
        P(i,j+2) =  patch[qq]->Uq*patch[qq]->cosalpha;
        P(i,j+3) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        Q(i,j)   =  patch[qq]->Uq*patch[qq]->sinalpha;
        Q(i,j+1) =  patch[qq]->Uq*patch[qq]->sinalpha;
        Q(i,j+2) =  patch[qq]->Uq*patch[qq]->sinalpha;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        P(i,j-1) =  patch[qq]->Uq*patch[qq]->cosalpha;
        P(i,j-2) =  patch[qq]->Uq*patch[qq]->cosalpha;
        P(i,j-3) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        Q(i,j-1) =  patch[qq]->Uq*patch[qq]->sinalpha;
        Q(i,j-2) =  patch[qq]->Uq*patch[qq]->sinalpha;
        Q(i,j-3) =  patch[qq]->Uq*patch[qq]->sinalpha;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        P(i,j)   =  patch[qq]->Uq*patch[qq]->sinalpha;
        P(i+1,j) =  patch[qq]->Uq*patch[qq]->sinalpha;
        P(i+2,j) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        Q(i+1,j) =  patch[qq]->Uq*patch[qq]->cosalpha;
        Q(i+2,j) =  patch[qq]->Uq*patch[qq]->cosalpha;
        Q(i+3,j) =  patch[qq]->Uq*patch[qq]->cosalpha;
        }
    }
    
    
    
     // Velocity components
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->velcomp_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        P(i-1,j) =  patch[qq]->U;
        P(i-2,j) =  patch[qq]->U;
        P(i-3,j) =  patch[qq]->U;
        
        Q(i,j)   =  patch[qq]->V;
        Q(i-1,j) =  patch[qq]->V;
        Q(i-2,j) =  patch[qq]->V;
        Q(i-3,j) =  patch[qq]->V;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        P(i,j)   =  patch[qq]->U;
        P(i,j+1) =  patch[qq]->U;
        P(i,j+2) =  patch[qq]->U;
        P(i,j+3) =  patch[qq]->U;
        
        Q(i,j)   =  patch[qq]->V;
        Q(i,j+1) =  patch[qq]->V;
        Q(i,j+2) =  patch[qq]->V;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        P(i,j)   =  patch[qq]->U;
        P(i,j-1) =  patch[qq]->U;
        P(i,j-2) =  patch[qq]->U;
        P(i,j-3) =  patch[qq]->U;
        
        Q(i,j-1) =  patch[qq]->V;
        Q(i,j-2) =  patch[qq]->V;
        Q(i,j-3) =  patch[qq]->V;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        P(i,j)   =  patch[qq]->U;
        P(i+1,j) =  patch[qq]->U;
        P(i+2,j) =  patch[qq]->U;
        
        Q(i,j)   =  patch[qq]->V;
        Q(i+1,j) =  patch[qq]->V;
        Q(i+2,j) =  patch[qq]->V;
        Q(i+3,j) =  patch[qq]->V;
        }
    
    }
    
}

void patchBC_2D::patchBC_rkioflow2D(lexer *p, ghostcell *pgc, slice &P, slice &Q, slice &U, slice &V)
{
    // Velocity components
    for(qq=0;qq<obj_count;++qq)
    //if(patch[qq]->velcomp_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        P(i-1,j) =  U(i-1,j);
        P(i-2,j) =  U(i-2,j);
        P(i-3,j) =  U(i-3,j);
        
        Q(i-1,j) =  V(i-1,j);
        Q(i-2,j) =  V(i-2,j);
        Q(i-3,j) =  V(i-3,j);
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        P(i,j+1) =  U(i,j+1);
        P(i,j+2) =  U(i,j+2);
        P(i,j+3) =  U(i,j+3);
        
        Q(i,j)   =  V(i,j);
        Q(i,j+1) =  V(i,j+1);
        Q(i,j+2) =  V(i,j+2);
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        P(i,j-1) =  U(i,j-1);
        P(i,j-2) =  U(i,j-2);
        P(i,j-3) =  U(i,j-3);
        
        Q(i,j-1) =  V(i,j-1);
        Q(i,j-2) =  V(i,j-2);
        Q(i,j-3) =  V(i,j-3);
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        P(i,j)   =  U(i,j);
        P(i+1,j) =  U(i,j+1);
        P(i+2,j) =  U(i,j+2);
        
        Q(i+1,j) =  V(i+1,j);
        Q(i+2,j) =  V(i+2,j);
        Q(i+3,j) =  V(i+3,j);
        }
    }
}


void patchBC_2D::patchBC_ioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    
    


} 

void patchBC_2D::patchBC_rkioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    
}



