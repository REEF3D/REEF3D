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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_f::wetdry(lexer* p, fdm_nhf* d, ghostcell* pgc, double *UH, double *VH, double *WH, slice &WL)
{
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    //------------
     
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
        if(WL(i,j)<=p->A544)
        {
        temp[IJ]=0;
        d->eta(i,j) = p->A544 - d->depth(i,j);
        WL(i,j) = d->eta(i,j) + d->depth(i,j);
        //cout<<d->eta(i,j)<<" "<<WL(i,j)<<endl;
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];
    
    //------------
    /*
    SLICELOOP4
    {
        if(WL(i,j)>=p->A544)
        p->wet[IJ]=1;
               
        if(WL(i,j)<p->A544)
        {
        p->wet[IJ]=0;
        d->eta(i,j) = p->A544 - d->depth(i,j);
        WL(i,j) = d->eta(i,j) + d->depth(i,j);
        }
    }*/
    
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
        
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    wetdry_fluxes(p,d,pgc,WL,d->U,d->V,d->W,UH,VH,WH);
    
    //
    SLICELOOP4
    p->deep[IJ]=p->wet[IJ];
    
    SLICELOOP4
    {
    if(p->wet[Ip1J]==0 || p->wet[Ip2J]==0 || p->wet[Ip3J]==0 || p->wet[Im1J]==0 || p->wet[Im2J]==0 || p->wet[Im3J]==0)
    p->deep[IJ]=0;
    
    if(p->wet[IJp1]==0 || p->wet[IJp2]==0 || p->wet[IJp3]==0 || p->wet[IJm1]==0 || p->wet[IJm2]==0 || p->wet[IJm3]==0)
    p->deep[IJ]=0;
    
    if(p->wet[Ip1Jp1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Im1Jm1]==0)
    p->deep[IJ]=0;
    }
    
    SLICELOOP4
    if(WL(i,j)<=10.0*p->A544)
    p->deep[IJ]=0;

    pgc->gcsl_start4Vint(p,p->deep,50);
    
}

void nhflow_fsf_f::wetdry_fluxes(lexer* p, fdm_nhf* d, ghostcell* pgc, slice &WL, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    
    // eta + WL
    SLICELOOP1  
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
        //d->ETAs(i,j) = d->eta(i,j);
        d->ETAn(i,j) = d->eta(i,j);
        //d->Ds(i,j) = WL(i,j);
        d->Dn(i,j) = WL(i,j);
        d->dfx(i,j) = d->depth(i,j);
        //d->dfx(i+1,j) = d->depth(i+1,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[Ip1J]==1)
        {
        d->ETAs(i,j) = d->eta(i+1,j);
        //d->ETAn(i,j) = d->eta(i+1,j);
        d->Ds(i,j) = WL(i+1,j);
        //d->Dn(i,j) = WL(i+1,j);
        d->dfx(i,j) = d->depth(i+1,j);
        //d->dfx(i-1,j) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[Ip1J]==0) 
        {
        d->ETAs(i,j) = d->eta(i,j);
        d->Ds(i,j) = WL(i,j);
        
        d->ETAn(i,j) = d->eta(i+1,j);
        d->Dn(i,j) = WL(i+1,j);
        }
        
        
        /*
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
        d->ETAs(i,j) = d->eta(i,j);
        d->ETAn(i,j) = d->eta(i,j);
        d->Ds(i,j) = WL(i,j);
        d->Dn(i,j) = WL(i,j);
        d->Dn(i-1,j) = WL(i,j);
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
        }*/
    }
    
    if(p->j_dir==1)
    SLICELOOP2 
    {
        if(p->wet[IJ]==1 && p->wet[IJp1]==0)
        {
        //d->ETAe(i,j) = d->eta(i,j);
        d->ETAw(i,j) = d->eta(i,j);
        //d->De(i,j) = WL(i,j);
        d->Dw(i,j) = WL(i,j);
        d->dfy(i,j) = d->depth(i,j);
        //d->dfy(i,j+1) = d->depth(i,j+1);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==1)
        {
        d->ETAe(i,j) = d->eta(i,j+1);
        //d->ETAw(i,j) = d->eta(i,j+1);
        d->De(i,j) = WL(i,j+1);
        //d->Dw(i,j) = WL(i,j+1);
        d->dfy(i,j) = d->depth(i,j+1);
        //d->dfy(i,j-1) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==0)
        {
        d->ETAe(i,j) = d->eta(i,j);
        d->De(i,j) = WL(i,j);
        
        d->ETAw(i,j) = d->eta(i,j+1);
        d->Dw(i,j) = WL(i,j+1);
        }
        
        /*
        if(p->wet[IJ]==1 && p->wet[IJp1]==0)
        {
        d->ETAe(i,j) = d->eta(i,j);
        d->ETAw(i,j) = d->eta(i,j+1);
        d->De(i,j) = WL(i,j);
        d->Dw(i,j) = WL(i,j+1);
        d->dfy(i,j) = d->depth(i,j);
        }
        
        else
        if(p->wet[IJ]==0 && p->wet[IJp1]==1)
        {
        d->ETAe(i,j) = d->eta(i,j);
        d->ETAw(i,j) = d->eta(i,j+1);
        d->De(i,j) = WL(i,j);
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
        }*/
    }
    
    // U,UH
    ULOOP 
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
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
        if(p->wet[IJ]==0 && p->wet[Ip1J]==1)
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
        if(p->wet[IJ]==0 && p->wet[Ip1J]==0)
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
        if(p->wet[IJ]==1 && p->wet[IJp1]==0)
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
        if(p->wet[IJ]==0 && p->wet[IJp1]==1)
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
        if(p->wet[IJ]==0 && p->wet[IJp1]==0)
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