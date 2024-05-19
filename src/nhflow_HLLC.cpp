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

#include"nhflow_HLLC.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"
#include"nhflow_flux_build_f.h"

nhflow_HLLC::nhflow_HLLC (lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_build_f(p,pgc,pBC);
}

nhflow_HLLC::~nhflow_HLLC()
{
}

void nhflow_HLLC::precalc(lexer* p, fdm_nhf* d, int ipol, slice &eta)
{
}

void nhflow_HLLC::start(lexer *&p, fdm_nhf *&d, int ipol, slice &eta)
{
    if(ipol==1)
    aij_U(p,d,1);

    if(ipol==2 && p->j_dir==1)
    aij_V(p,d,2);

    if(ipol==3)
    aij_W(p,d,3);
        
    if(ipol==4)
    aij_E(p,d,4);
}

void nhflow_HLLC::aij_U(lexer* p,fdm_nhf* d, int ipol)
{
    // HLLC flux 
    pflux->start_U(p,d,pgc);
    HLLC(p,d,d->UHs,d->UHn,d->UHe,d->UHw,d->SSx,d->SSx,d->Ue,d->Uw);
    
    pgc->start1V(p,d->Fx,10);
    pgc->start2V(p,d->Fy,10);
    pgc->start3V(p,d->Fz,10);
    
    LOOP
    WETDRY
    {
    d->F[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

void nhflow_HLLC::aij_V(lexer* p, fdm_nhf* d, int ipol)
{
    // HLLC flux 
    if(p->j_dir==1)
    {
    pflux->start_V(p,d,pgc);
    HLLC(p,d,d->VHs,d->VHn,d->VHe,d->VHw,d->Vs,d->Vn,d->SSy,d->SSy);
    
    LOOP
    WETDRY
    {
    if(p->wet[IJp1]==0 && p->flag2[IJp1K]>0)
    d->Fy[IJK] = 0.0;
    
    if(p->wet[IJm1]==0 && p->flag2[IJm1K]>0)
    d->Fy[IJm1K] = 0.0;
    }
    
    pgc->start1V(p,d->Fx,11);
    pgc->start2V(p,d->Fy,11);
    pgc->start3V(p,d->Fz,11);
    
    LOOP
    WETDRY
    {
    d->G[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }  

    }  
}

void nhflow_HLLC::aij_W(lexer* p,fdm_nhf* d, int ipol)
{
    // HLLC flux 
    pflux->start_W(p,d,pgc);
    HLLC(p,d,d->WHs,d->WHn,d->WHe,d->WHw,d->Ws,d->Wn,d->We,d->Ww);
    //HLL(p,d,d->WHs,d->WHn,d->WHe,d->WHw);
    
    pgc->start1V(p,d->Fx,12);
    pgc->start2V(p,d->Fy,12);
    pgc->start3V(p,d->Fz,12);
    
    LOOP
    WETDRY
    {
    d->H[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

void nhflow_HLLC::aij_E(lexer* p,fdm_nhf* d, int ipol)
{
    // HLLC flux 
    pflux->start_E(p,d,pgc);
    HLLC_E(p,d);
    
    LOOP
    WETDRY
    {
    if(p->wet[Ip1J]==0)
    d->Fx[IJK] = 0.0;
    
    if(p->wet[Im1J]==0)
    d->Fx[Im1JK] = 0.0;
    
    if(p->wet[IJp1]==0)
    d->Fy[IJK] = 0.0;
    
    if(p->wet[IJm1]==0)
    d->Fy[IJm1K] = 0.0;
    }
    
    pgc->start1V(p,d->Fx,14);
    pgc->start2V(p,d->Fy,14); 
}

void nhflow_HLLC::HLLC(lexer* p,fdm_nhf* d, double *Us, double *Un, double *Ue, double *Uw, double *SSxs, double *SSxn, double *SSye, double *SSyw)
{
    // HLLC flux
    ULOOP
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==1 && p->wet[Im1J]==1 && p->wet[Ip2J]==1)
        {
            FsS = d->Ds(i,j)*(d->Ss[IJK] - d->Us[IJK] + 1.0e-10)/(d->Ss[IJK] - d->SSx[IJK] + 1.0e-10)*SSxs[IJK];
            FnS = d->Dn(i,j)*(d->Sn[IJK] - d->Un[IJK] + 1.0e-10)/(d->Sn[IJK] - d->SSx[IJK] + 1.0e-10)*SSxn[IJK];
     
            if(d->Ss[IJK]>=0.0)
            d->Fx[IJK] = d->Fs[IJK];
            
            else
            if(d->Sn[IJK]<=0.0)
            d->Fx[IJK] = d->Fn[IJK];
            
            else
            if(d->SSx[IJK]>=0.0)
            d->Fx[IJK] = d->Fs[IJK] + d->Ss[IJK]*(FsS - Us[IJK]);
            
            else
            d->Fx[IJK] = d->Fn[IJK] + d->Sn[IJK]*(FnS - Un[IJK]);
        }
        
        if(p->wet[IJ]==0 || p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[Ip2J]==0)
        {
            if(d->Ss[IJK]>=0.0)
            d->Fx[IJK] = d->Fs[IJK];
            
            else
            if(d->Sn[IJK]<=0.0)
            d->Fx[IJK] = d->Fn[IJK];
            
            else
            {
            denom = d->Sn[IJK]-d->Ss[IJK];
            denom = fabs(denom)>1.0e-10?denom:1.0e10;
            
            d->Fx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(Un[IJK] - Us[IJK]))/denom;
            }
        }
    }
        
    // HLLC flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
        if(p->wet[IJ]==1 && p->wet[IJp1]==1 && p->wet[IJm1]==1 && p->wet[IJp2]==1)
        {
            FeS = d->De(i,j)*(d->Se[IJK] - d->Ve[IJK] + 1.0e-10)/(d->Se[IJK] - d->SSy[IJK] + 1.0e-10)*SSye[IJK];
            FwS = d->Dw(i,j)*(d->Sw[IJK] - d->Vw[IJK] + 1.0e-10)/(d->Sw[IJK] - d->SSy[IJK] + 1.0e-10)*SSyw[IJK];
     
            if(d->Se[IJK]>=0.0)
            d->Fy[IJK] = d->Fe[IJK];
            
            else
            if(d->Sw[IJK]<=0.0)
            d->Fy[IJK] = d->Fw[IJK];
            
            else
            if(d->SSy[IJK]>=0.0)
            d->Fy[IJK] = d->Fe[IJK] + d->Se[IJK]*(FeS - Ue[IJK]);
            
            else
            d->Fy[IJK] = d->Fw[IJK] + d->Sw[IJK]*(FwS - Uw[IJK]);
        }
        
        if(p->wet[IJ]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0 || p->wet[IJp2]==0)
        {
            if(d->Se[IJK]>=0.0)
            d->Fy[IJK] = d->Fe[IJK];
            
            else
            if(d->Sw[IJK]<=0.0)
            d->Fy[IJK] = d->Fw[IJK];
            
            else
            {
            denom = d->Sw[IJK]-d->Se[IJK];
            denom = fabs(denom)>1.0e-10?denom:1.0e10;
            
            d->Fy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(Uw[IJK] - Ue[IJK]))/denom;
            }
        }


    }
}

void nhflow_HLLC::HLLC_E(lexer* p,fdm_nhf* d)
{
    
    // HLLC flux
    ULOOP
    {
        if(p->wet[IJ]==1 && p->wet[Ip1J]==1 && p->wet[Im1J]==1 && p->wet[Ip2J]==1)
        {
            
            FsS = d->Ds(i,j)*(d->Ss[IJK] - d->Us[IJK] + 1.0e-10)/(d->Ss[IJK] - d->SSx[IJK] + 1.0e-10);
            FnS = d->Dn(i,j)*(d->Sn[IJK] - d->Un[IJK] + 1.0e-10)/(d->Sn[IJK] - d->SSx[IJK] + 1.0e-10);
     
            if(d->Ss[IJK]>=0.0)
            d->Fx[IJK] = d->Fs[IJK];
            
            else
            if(d->Sn[IJK]<=0.0)
            d->Fx[IJK] = d->Fn[IJK];
            
            else
            if(d->SSx[IJK]>=0.0)
            d->Fx[IJK] = d->Fs[IJK] + d->Ss[IJK]*(FsS - d->Ds(i,j));
            
            else
            d->Fx[IJK] = d->Fn[IJK] + d->Sn[IJK]*(FnS - d->Dn(i,j));
        }
        
        if(p->wet[IJ]==0 || p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[Ip2J]==0)
        {
            
            if(d->Ss[IJK]>=0.0)
            d->Fx[IJK] = d->Fs[IJK];
            
            else
            if(d->Sn[IJK]<=0.0)
            d->Fx[IJK] = d->Fn[IJK];
            
            else
            {
            denom = d->Sn[IJK]-d->Ss[IJK];
            denom = fabs(denom)>1.0e-10?denom:1.0e10;
            
            d->Fx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(d->Dn(i,j) - d->Ds(i,j)))/denom;
            }
        }
    }
    
    // HLLC flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
        if(p->wet[IJ]==1 && p->wet[IJp1]==1 && p->wet[IJm1]==1 && p->wet[IJp2]==1)
        {
            FeS = d->De(i,j)*(d->Se[IJK] - d->Ve[IJK] + 1.0e-10)/(d->Se[IJK] - d->SSy[IJK] + 1.0e-10);
            FwS = d->Dw(i,j)*(d->Sw[IJK] - d->Vw[IJK] + 1.0e-10)/(d->Sw[IJK] - d->SSy[IJK] + 1.0e-10);
     
            if(d->Se[IJK]>=0.0)
            d->Fy[IJK] = d->Fe[IJK];
            
            else
            if(d->Sw[IJK]<=0.0)
            d->Fy[IJK] = d->Fw[IJK];
            
            else
            if(d->SSy[IJK]>=0.0)
            d->Fy[IJK] = d->Fe[IJK] + d->Se[IJK]*(FeS - d->De(i,j));
            
            else
            d->Fy[IJK] = d->Fw[IJK] + d->Sw[IJK]*(FwS - d->Dw(i,j));
        }
        
        if(p->wet[IJ]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0 || p->wet[IJp2]==0)
        {
            if(d->Se[IJK]>=0.0)
            d->Fy[IJK] = d->Fe[IJK];
            
            else
            if(d->Sw[IJK]<=0.0)
            d->Fy[IJK] = d->Fw[IJK];
            
            else
            {
            denom = d->Sw[IJK]-d->Se[IJK];
            denom = fabs(denom)>1.0e-10?denom:1.0e10;
            
            d->Fy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(d->Dw(i,j) - d->De(i,j)))/denom;
            }
        }
    }
}

void nhflow_HLLC::HLL(lexer *p,fdm_nhf *&d, double *Us, double *Un, double *Ue, double *Uw)
{
    // HLL flux
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        d->Fx[IJK] = d->Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        d->Fx[IJK] = d->Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->Fx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(Un[IJK] - Us[IJK]))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    {
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        d->Fy[IJK] = d->Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        d->Fy[IJK] = d->Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->Fy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(Uw[IJK] - Ue[IJK]))/denom;
        }
    }
    }
}
