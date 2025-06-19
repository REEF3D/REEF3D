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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_f::wetdry(lexer* p, fdm_nhf* d, ghostcell* pgc, double *UH, double *VH, double *WH, slice &WL)
{
        if(p->count==0)
        {
            SLICELOOP4
            if(WL(i,j)<=p->A544+eps)
            {
            temp[IJ]=0;
            d->eta(i,j) =  p->A544 - d->depth(i,j);
            WL(i,j) = p->A544;
            }
            
            
            SLICEBASELOOP
            if(p->flagslice4[IJ]<0)
            {
            p->wet[IJ]=0;
            temp[IJ]=0;
            }
        pgc->gcsl_start4(p,d->eta,gcval_eta);
        }
    
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    // ----------------
    // ----------------
    if(p->A540==1)
    {
    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
            if(p->wet[Ip1J]==1 && d->eta(i,j)<d->eta(i+1,j) && WL(i+1,j)>p->A544+eps)
            temp[IJ]=1;
            
            if(p->wet[Im1J]==1 && d->eta(i,j)<d->eta(i-1,j) && WL(i-1,j)>p->A544+eps)
            temp[IJ]=1;
            
            if(p->wet[IJp1]==1 && d->eta(i,j)<d->eta(i,j+1) && WL(i,j+1)>p->A544+eps && p->j_dir==1)
            temp[IJ]=1;
            
            if(p->wet[IJm1]==1 && d->eta(i,j)<d->eta(i,j-1) && WL(i,j-1)>p->A544+eps && p->j_dir==1)
            temp[IJ]=1;
        }
        
        else              
        if(WL(i,j)<=p->A544+eps)
        {
        temp[IJ]=0;
        d->eta(i,j) =  p->A544 - d->depth(i,j);
        WL(i,j) = p->A544;
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];
     }
     
     pgc->gcsl_start4(p,d->eta,gcval_eta);
     pgc->gcsl_start4Vint(p,p->wet,50);
    
    // ----------------
    // ----------------
    if(p->A540==2)
    { 
    
        SLICELOOP4
        {
            if(WL(i,j)>=p->A544)
            p->wet[IJ]=1;
            
            if(WL(i,j)<=p->A544+eps)
            {
            p->wet[IJ]=0;
            d->eta(i,j) =  p->A544 - d->depth(i,j);
            WL(i,j) = p->A544;
            }
        }
    }
    
    // avoid isolated wetdry
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
    if(p->wet[Im1J]==0 && p->flagslice4[Ip1J]<0)
    p->wet[IJ]=0;
    
    if(p->wet[Ip1J]==0 && p->flagslice4[Im1J]<0)
    p->wet[IJ]=0;
    
    if(p->wet[IJm1]==0 && p->flagslice4[IJp1]<0 && p->j_dir==1)
    p->wet[IJ]=0;
    
    if(p->wet[IJp1]==0 && p->flagslice4[IJm1]<0 && p->j_dir==1)
    p->wet[IJ]=0;
    
    
    
    if(p->wet[Im1J]==0 && p->wet[Ip1J]==0)
    p->wet[IJ]=0;
    
    if(p->wet[Ip1J]==0 && p->wet[Im1J]==0)
    p->wet[IJ]=0;
    
    if(p->wet[IJm1]==0 && p->wet[IJp1]==0 && p->j_dir==1)
    p->wet[IJ]=0;
    
    if(p->wet[IJp1]==0 && p->wet[IJm1]==0 && p->j_dir==1)
    p->wet[IJ]=0;
        
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    //------------
    LOOP
    if(p->wet[IJ]==0)
    {
        d->U[IJK] = 0.0;
        d->V[IJK] = 0.0;
        d->W[IJK] = 0.0;
        
        UH[IJK] = 0.0;
        VH[IJK] = 0.0;
        WH[IJK] = 0.0;
    }
    
    FLOOP
    if(p->wet[IJ]==0)
    d->omegaF[FIJK] = 0.0;
    
        
    // --------------------------
    wetdry_fluxes(p,d,pgc,WL,d->U,d->V,d->W,UH,VH,WH);
    // --------------------------
    
    //
    SLICELOOP4
    p->deep[IJ]=p->wet[IJ];
    
    SLICELOOP4
    {
    if(p->wet[Ip1J]==0 || p->wet[Ip2J]==0 || p->wet[Im1J]==0 || p->wet[Im2J]==0)
    p->deep[IJ]=0;
    
    
    if(p->j_dir==1)
    if(p->wet[IJp1]==0 || p->wet[IJp2]==0 || p->wet[IJp3]==0 || p->wet[IJm1]==0 || p->wet[IJm2]==0 || p->wet[IJm3]==0)
    p->deep[IJ]=0;
    
    if(p->j_dir==1)
    if(p->wet[Ip1Jp1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Im1Jm1]==0)
    p->deep[IJ]=0;
    }
    
    SLICELOOP4
    if(WL(i,j)<=p->A545*p->A544)
    p->deep[IJ]=0;

    pgc->gcsl_start4Vint(p,p->deep,50);
}