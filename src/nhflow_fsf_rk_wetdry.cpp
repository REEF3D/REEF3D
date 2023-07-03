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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_rk::wetdry(lexer* p, fdm_nhf* d, ghostcell* pgc, double *UH, double *VH, double *WH, slice &WL)
{
    
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
     
    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
            if(p->wet[Ip1J]==1 && d->eta(i,j)<d->eta(i+1,j) && d->WL(i+1,j)>p->A544+eps)
            temp[IJ]=1;
            
            if(p->wet[Im1J]==1 && d->eta(i,j)<d->eta(i-1,j) && d->WL(i-1,j)>p->A544+eps)
            temp[IJ]=1;
            
            if(p->wet[IJp1]==1 && d->eta(i,j)<d->eta(i,j+1) && d->WL(i,j+1)>p->A544+eps && p->j_dir==1)
            temp[IJ]=1;
            
            if(p->wet[IJm1]==1 && d->eta(i,j)<d->eta(i,j-1) && d->WL(i,j-1)>p->A544+eps && p->j_dir==1)
            temp[IJ]=1;
        }
        
        else              
        if(d->WL(i,j)<=p->A544)
        {
        temp[IJ]=0;
        d->eta(i,j) = p->A544 - d->depth(i,j) - eps;
        WL(i,j) = d->eta(i,j) + d->depth(i,j) + eps;
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];

    
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
    
    /*
    SLICELOOP4
    if(p->wet[IJ]==0)
    //if(d->eta(i,j)< -p->wd  + d->bed(i,j) - p->A544 + 1.0e-20)
    {
    //d->eta(i,j) = -p->wd  + d->bed(i,j) - p->A544 - eps;
    
    d->eta(i,j) = p->A544 - d->depth(i,j) - eps;
    d->WL(i,j) = d->eta(i,j) + d->depth(i,j) + eps;
    }*/
    
        
    pgc->gcsl_start4(p,d->eta,1);
    pgc->gcsl_start4Vint(p,p->wet,50);
}


void nhflow_fsf_rk::wetdry_fluxes(lexer* p, fdm_nhf* d, ghostcell* pgc, slice &WL, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    // eta + WL
    SLICELOOP1  
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
        d->ETAs(i,j) = d->eta(i,j);

        d->Ds(i,j) = WL(i,j);
        
        d->dfx(i,j) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==1 && p->wet[Im1J]==0)
        {
        d->ETAn(i-1,j) = d->eta(i,j);

        d->Dn(i-1,j) = WL(i,j);
        
        d->dfx(i-1,j) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0)
        {
        d->ETAs(i,j) = d->eta(i,j);

        d->Ds(i,j) = WL(i,j);
        
        
        d->ETAn(i-1,j) = d->eta(i,j);

        d->Dn(i-1,j) = WL(i,j);
        }
    }
    
    if(p->j_dir==1)
    SLICELOOP2 
    {
        if(p->wet[IJ]==1 && p->wet[IJp1]==0)
        {
        d->ETAe(i,j) = -p->wd  + d->bed(i,j) - p->A544 - 1.0e-15; 
        
        d->De(i,j) = d->WL(i,j);
        }
        
        else
        if(p->wet[IJ]==1 && p->wet[IJm1]==0)
        {
        d->ETAw(i,j-1) = -p->wd  + d->bed(i,j) - p->A544 - 1.0e-15;
        }
        
        else
        if(p->wet[IJ]==0)
        {
        d->ETAe(i,j) = -p->wd  + d->bed(i,j) - p->A544 - 1.0e-15;
        d->ETAw(i,j-1) = -p->wd  + d->bed(i,j) - p->A544 - 1.0e-15;
        }
    }
    

    // U,UH
    // reconstruct
    ULOOP 
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
        d->Us[IJK] = 0.0;
        d->Vs[IJK] = 0.0;
        d->Ws[IJK] = 0.0;
        
        d->UHs[IJK] = 0.0;
        d->VHs[IJK] = 0.0;
        d->WHs[IJK] = 0.0;
        }
        
        else
        if(p->wet[IJ]==1 && p->wet[Im1J]==0)
        {
        d->Un[Im1JK] = 0.0;
        d->Vn[Im1JK] = 0.0;
        d->Wn[Im1JK] = 0.0;
        
        d->UHn[Im1JK] = 0.0;
        d->VHn[Im1JK] = 0.0;
        d->WHn[Im1JK] = 0.0;
        }
        
        else
        if(p->wet[IJ]==0)
        {
        d->Us[IJK] = 0.0;
        d->Vs[IJK] = 0.0;
        d->Ws[IJK] = 0.0;
        
        d->UHs[IJK] = 0.0;
        d->VHs[IJK] = 0.0;
        d->WHs[IJK] = 0.0;
        
        d->Un[Im1JK] = 0.0;
        d->Vn[Im1JK] = 0.0;
        d->Wn[Im1JK] = 0.0;
        
        d->UHn[Im1JK] = 0.0;
        d->VHn[Im1JK] = 0.0;
        d->WHn[Im1JK] = 0.0;
        }
    }
    
    
}