/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
        cout<<"Ui_patch_0: "<<u(i,j,k)<<" "<<u(i+1,j,k)<<" "<<u(i+2,j,k)<<endl;
        u(i,j,k)   =  patch[qq]->Uio;
        u(i+1,j,k) =  patch[qq]->Uio;
        u(i+2,j,k) =  patch[qq]->Uio;
        
        v(i+1,j,k) =  0.0;
        v(i+2,j,k) =  0.0;
        v(i+3,j,k) =  0.0;
        
        w(i+1,j,k) =  0.0;
        w(i+2,j,k) =  0.0;
        w(i+3,j,k) =  0.0;
        cout<<"Ui_patch_1: "<<u(i,j,k)<<" "<<u(i+1,j,k)<<" "<<u(i+2,j,k)<<endl;
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
        cout<<"Uq_patch_0: "<<u(i,j,k)<<" "<<u(i+1,j,k)<<" "<<u(i+2,j,k)<<endl;
        u(i,j,k)   =  patch[qq]->Uq*patch[qq]->sinalpha;
        u(i+1,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        u(i+2,j,k) =  patch[qq]->Uq*patch[qq]->sinalpha;
        
        v(i+1,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        v(i+2,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        v(i+3,j,k) =  patch[qq]->Uq*patch[qq]->cosalpha;
        
        w(i+1,j,k) =  0.0;
        w(i+2,j,k) =  0.0;
        w(i+3,j,k) =  0.0;
        
        cout<<"Uq_patch_1: "<<u(i,j,k)<<" "<<u(i+1,j,k)<<" "<<u(i+2,j,k)<<endl;
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
        cout<<"Ucomp_patch_0: "<<u(i-1,j,k)<<" "<<u(i,j,k)<<" "<<u(i+1,j,k)<<" "<<u(i+2,j,k)<<endl;
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
        cout<<"Ucomp_patch_1: "<<u(i-1,j,k)<<" "<<u(i,j,k)<<" "<<u(i+1,j,k)<<" "<<u(i+2,j,k)<<endl;
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


void patchBC::patchBC_ioflow2D(lexer *p, ghostcell*, slice&, slice&, slice&, slice&)
{
    cout<<p->mpirank<<" patchBC_2D: inflow 3D "<<endl;
}


