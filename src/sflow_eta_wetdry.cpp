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

#include"sflow_eta.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sflow_eta::wetdry(lexer* p, fdm2D* b, ghostcell* pgc, slice &WL)
{
    // A 243  0: no wetting-drying, 1: depth criterion, 2: depth criterion + wetting from higher wet neighbours
    
    if(p->count==0)
    {
        SLICELOOP4
        if(WL(i,j)<=wd_criterion+eps && p->A243>0)
        {
        temp[IJ]=0;
        b->eta(i,j) = wd_criterion - b->depth(i,j);
        WL(i,j) = wd_criterion;
        }
        
        SLICEBASELOOP
        if(p->flagslice4[IJ]<0)
        {
        p->wet[IJ]=0;
        temp[IJ]=0;
        }
    }
    
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    // ----------------
    if(p->A243==0)
    {
        SLICELOOP4
        {
        temp[IJ]=1;
        
        if(WL(i,j)<wd_criterion)
        {
        b->eta(i,j) = wd_criterion - b->depth(i,j);
        WL(i,j) = wd_criterion;
        }
        }
    }
    
    if(p->A243==1)
    {
        SLICELOOP4
        {
            if(WL(i,j)>=wd_criterion+eps)
            temp[IJ]=1;
            
            if(WL(i,j)<wd_criterion+eps)
            {
            temp[IJ]=0;
            b->eta(i,j) = wd_criterion - b->depth(i,j);
            WL(i,j) = wd_criterion;
            }
        }
    }
    
    if(p->A243==2)
    {
        SLICELOOP4
        {
            if(p->wet[IJ]==0)
            {
                if(p->wet[Ip1J]==1 && b->eta(i,j)<b->eta(i+1,j) && WL(i+1,j)>wd_criterion+eps)
                temp[IJ]=1;
                
                if(p->wet[Im1J]==1 && b->eta(i,j)<b->eta(i-1,j) && WL(i-1,j)>wd_criterion+eps)
                temp[IJ]=1;
                
                if(p->wet[IJp1]==1 && b->eta(i,j)<b->eta(i,j+1) && WL(i,j+1)>wd_criterion+eps && p->j_dir==1)
                temp[IJ]=1;
                
                if(p->wet[IJm1]==1 && b->eta(i,j)<b->eta(i,j-1) && WL(i,j-1)>wd_criterion+eps && p->j_dir==1)
                temp[IJ]=1;
            }
            
            else              
            if(WL(i,j)<=wd_criterion+eps)
            {
            temp[IJ]=0;
            b->eta(i,j) = wd_criterion - b->depth(i,j);
            WL(i,j) = wd_criterion;
            }
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    // avoid isolated wet cells
    if(p->A243>0)
    {
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
        
        if(p->wet[IJm1]==0 && p->wet[IJp1]==0 && p->j_dir==1)
        p->wet[IJ]=0;
        }
        
        pgc->gcsl_start4Vint(p,p->wet,50);
    }
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    pgc->gcsl_start4(p,WL,gcval_eta);
    
    // dry cells carry no momentum
    SLICELOOP4
    if(p->wet[IJ]==0)
    {
    b->U(i,j) = 0.0;
    b->V(i,j) = 0.0;
    b->W(i,j) = 0.0;
    b->UH(i,j) = 0.0;
    b->VH(i,j) = 0.0;
    b->WH(i,j) = 0.0;
    }
    
    // face wet flags (used by the turbulence model)
    SLICELOOP1
    b->wet1(i,j) = (p->wet[IJ]==1 && p->wet[Ip1J]==1)?1:0;
    
    SLICELOOP2
    b->wet2(i,j) = (p->wet[IJ]==1 && p->wet[IJp1]==1)?1:0;
    
    pgc->gcsl_start1int(p,b->wet1,50);
    pgc->gcsl_start2int(p,b->wet2,50);
    
    // water depth for the other modules
    SLICELOOP4
    b->hp(i,j) = WL(i,j);
    
    pgc->gcsl_start4(p,b->hp,gcval_eta);
    
    // deep flag: non-hydrostatic pressure only away from the shoreline (A 221 0)
    wetdrydeep(p,b,pgc,WL);
    
    // gcslin update
    if(p->count<=1)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    if(p->wet[IJ]==0)
    p->gcslin[n][5]=0;
    }
}

void sflow_eta::wetdry_fluxes(lexer* p, fdm2D* b, ghostcell* pgc, slice &WL)
{
    if(p->A243==0)
    return;
    
    // x-faces
    SLICELOOP1
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
        b->ETAs(i,j) = b->eta(i,j);
        b->ETAn(i,j) = b->eta(i,j);
        b->Ds(i,j) = WL(i,j);
        b->Dn(i,j) = WL(i,j);
        b->dfx(i,j) = b->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[Ip1J]==1)
        {
        b->ETAs(i,j) = b->eta(i+1,j);
        b->ETAn(i,j) = b->eta(i+1,j);
        b->Ds(i,j) = WL(i+1,j);
        b->Dn(i,j) = WL(i+1,j);
        b->dfx(i,j) = b->depth(i+1,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[Ip1J]==0) 
        {
        b->ETAs(i,j) = b->eta(i,j);
        b->Ds(i,j) = WL(i,j);
        b->ETAn(i,j) = b->eta(i+1,j);
        b->Dn(i,j) = WL(i+1,j);
        }
        
        if(!(p->wet[IJ]==1 && p->wet[Ip1J]==1))
        {
        b->Us(i,j) = 0.0;
        b->Un(i,j) = 0.0;
        b->UHs(i,j) = 0.0;
        b->UHn(i,j) = 0.0;
        b->VHs(i,j) = 0.0;
        b->VHn(i,j) = 0.0;
        b->WHs(i,j) = 0.0;
        b->WHn(i,j) = 0.0;
        b->QUs(i,j) = 0.0;
        b->QUn(i,j) = 0.0;
        b->QVs(i,j) = 0.0;
        b->QVn(i,j) = 0.0;
        }
    }
    
    // y-faces
    if(p->j_dir==1)
    SLICELOOP2
    {
        if(p->wet[IJ]==1 && p->wet[IJp1]==0)
        {
        b->ETAe(i,j) = b->eta(i,j);
        b->ETAw(i,j) = b->eta(i,j);
        b->De(i,j) = WL(i,j);
        b->Dw(i,j) = WL(i,j);
        b->dfy(i,j) = b->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==1)
        {
        b->ETAe(i,j) = b->eta(i,j+1);
        b->ETAw(i,j) = b->eta(i,j+1);
        b->De(i,j) = WL(i,j+1);
        b->Dw(i,j) = WL(i,j+1);
        b->dfy(i,j) = b->depth(i,j+1);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==0)
        {
        b->ETAe(i,j) = b->eta(i,j);
        b->De(i,j) = WL(i,j);
        b->ETAw(i,j) = b->eta(i,j+1);
        b->Dw(i,j) = WL(i,j+1);
        }
        
        if(!(p->wet[IJ]==1 && p->wet[IJp1]==1))
        {
        b->Ve(i,j) = 0.0;
        b->Vw(i,j) = 0.0;
        b->UHe(i,j) = 0.0;
        b->UHw(i,j) = 0.0;
        b->VHe(i,j) = 0.0;
        b->VHw(i,j) = 0.0;
        b->WHe(i,j) = 0.0;
        b->WHw(i,j) = 0.0;
        b->QUe(i,j) = 0.0;
        b->QUw(i,j) = 0.0;
        b->QVe(i,j) = 0.0;
        b->QVw(i,j) = 0.0;
        }
    }
    
    // consistent face depth for the bed-slope source across partitions
    pgc->gcsl_start1(p,b->dfx,1);
    pgc->gcsl_start2(p,b->dfy,1);
}
