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

void nhflow_fsf_f::wetdry_fluxes(lexer* p, fdm_nhf* d, ghostcell* pgc, slice &WL, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
   if(p->A540==1)
   { 
    // eta + WL
    SLICELOOP1
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
        d->ETAs(i,j) = d->eta(i,j);
        d->ETAn(i,j) = d->eta(i,j);
        d->Ds(i,j) = WL(i,j);
        d->Dn(i,j) = WL(i,j);
        d->dfx(i,j) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[Ip1J]==1)
        {
        d->ETAs(i,j) = d->eta(i+1,j);
        d->ETAn(i,j) = d->eta(i+1,j);
        d->Ds(i,j) = WL(i+1,j);
        d->Dn(i,j) = WL(i+1,j);
        d->dfx(i,j) = d->depth(i+1,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[Ip1J]==0) 
        {
        d->ETAs(i,j) = d->eta(i,j);
        d->Ds(i,j) = WL(i,j);
        
        d->ETAn(i,j) = d->eta(i+1,j);
        d->Dn(i,j) = WL(i+1,j);
        }
    }
    
    if(p->j_dir==1)
    SLICELOOP2
    {
        if(p->wet[IJ]==1 && p->wet[IJp1]==0)
        {
        d->ETAe(i,j) = d->eta(i,j);
        d->ETAw(i,j) = d->eta(i,j);
        d->De(i,j) = WL(i,j);
        d->Dw(i,j) = WL(i,j);
        d->dfy(i,j) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==1)
        {
        d->ETAe(i,j) = d->eta(i,j+1);
        d->ETAw(i,j) = d->eta(i,j+1);
        d->De(i,j) = WL(i,j+1);
        d->Dw(i,j) = WL(i,j+1);
        d->dfy(i,j) = d->depth(i,j+1);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==0)
        {
        d->ETAe(i,j) = d->eta(i,j);
        d->De(i,j) = WL(i,j);
        
        d->ETAw(i,j) = d->eta(i,j+1);
        d->Dw(i,j) = WL(i,j+1);
        }
    }
    
    // U,UH
    ULOOP
    {
        if((p->wet[IJ]==1 && p->wet[Ip1J]==0) || (p->DF[IJK]==1 && p->DF[Ip1JK]<0))
        {
        d->Un[IJK] = 0.0;
        d->Vn[IJK] = 0.0;
        d->Wn[IJK] = 0.0;
    
        d->UHn[IJK] = 0.0;
        d->VHn[IJK] = 0.0;
        d->WHn[IJK] = 0.0;
        
        d->Us[IJK] = 0.0;
        d->Vs[IJK] = 0.0;
        d->Ws[IJK] = 0.0;
        
        d->UHs[IJK] = 0.0;
        d->VHs[IJK] = 0.0;
        d->WHs[IJK] = 0.0;
        }
        
        else
        if((p->wet[IJ]==0 && p->wet[Ip1J]==1)  || (p->DF[IJK]<0 && p->DF[Ip1JK]==1))
        {
        d->Un[IJK] = 0.0;
        d->Vn[IJK] = 0.0;
        d->Wn[IJK] = 0.0;
    
        d->UHn[IJK] = 0.0;
        d->VHn[IJK] = 0.0;
        d->WHn[IJK] = 0.0;
        
        d->Us[IJK] = 0.0;
        d->Vs[IJK] = 0.0;
        d->Ws[IJK] = 0.0;
        
        d->UHs[IJK] = 0.0;
        d->VHs[IJK] = 0.0;
        d->WHs[IJK] = 0.0;
        }

        else
        if((p->wet[IJ]==0 && p->wet[Ip1J]==0)  || (p->DF[IJK]<0 && p->DF[Ip1JK]<0))
        {
        d->Un[IJK] = 0.0;
        d->Vn[IJK] = 0.0;
        d->Wn[IJK] = 0.0;
        
        d->UHn[IJK] = 0.0;
        d->VHn[IJK] = 0.0;
        d->WHn[IJK] = 0.0;
        
        d->Us[IJK] = 0.0;
        d->Vs[IJK] = 0.0;
        d->Ws[IJK] = 0.0;
        
        d->UHs[IJK] = 0.0;
        d->VHs[IJK] = 0.0;
        d->WHs[IJK] = 0.0;
        }
    }
    
    VLOOP
    {

        if((p->wet[IJ]==1 && p->wet[IJp1]==0) || (p->DF[IJK]==1 && p->DF[IJp1K]<0))
        {
        d->Uw[IJK] = 0.0;
        d->Vw[IJK] = 0.0;
        d->Ww[IJK] = 0.0;
    
        d->UHw[IJK] = 0.0;
        d->VHw[IJK] = 0.0;
        d->WHw[IJK] = 0.0;
        
        d->Ue[IJK] = 0.0;
        d->Ve[IJK] = 0.0;
        d->We[IJK] = 0.0;
        
        d->UHe[IJK] = 0.0;
        d->VHe[IJK] = 0.0;
        d->WHe[IJK] = 0.0;
        }
        
        else
        if((p->wet[IJ]==0 && p->wet[IJp1]==1) || (p->DF[IJK]<0 && p->DF[IJp1K]==1))
        {
        d->Uw[IJK] = 0.0;
        d->Vw[IJK] = 0.0;
        d->Ww[IJK] = 0.0;
    
        d->UHw[IJK] = 0.0;
        d->VHw[IJK] = 0.0;
        d->WHw[IJK] = 0.0;
        
        d->Ue[IJK] = 0.0;
        d->Ve[IJK] = 0.0;
        d->We[IJK] = 0.0;
        
        d->UHe[IJK] = 0.0;
        d->VHe[IJK] = 0.0;
        d->WHe[IJK] = 0.0;
        }
        
        else
        if((p->wet[IJ]==0 && p->wet[IJp1]==0) || (p->DF[IJK]<0 && p->DF[IJp1K]<0))
        {
        d->Uw[IJK] = 0.0;
        d->Vw[IJK] = 0.0;
        d->Ww[IJK] = 0.0;
        
        d->UHw[IJK] = 0.0;
        d->VHw[IJK] = 0.0;
        d->WHw[IJK] = 0.0;
        
        d->Ue[IJK] = 0.0;
        d->Ve[IJK] = 0.0;
        d->We[IJK] = 0.0;
        
        d->UHe[IJK] = 0.0;
        d->VHe[IJK] = 0.0;
        d->WHe[IJK] = 0.0;
        }
    }
   }
    
}