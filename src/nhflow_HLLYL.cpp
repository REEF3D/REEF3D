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

#include"nhflow_HLLYL.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"
#include"nhflow_flux_build_f.h"
#include"vrans.h"

nhflow_HLLYL::nhflow_HLLYL(lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_build_f(p,pgc,pBC);
}

nhflow_HLLYL::~nhflow_HLLYL()
{
}

void nhflow_HLLYL::precalc(lexer* p, fdm_nhf* d, int ipolL, slice &eta)
{
}

void nhflow_HLLYL::start(lexer *&p, fdm_nhf *&d, int ipol, slice &eta)
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

void nhflow_HLLYL::aij_U(lexer *&p,fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_U(p,d,pgc);
    HLL(p,d,d->UHs,d->UHn,d->UHe,d->UHw);
    
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

void nhflow_HLLYL::aij_V(lexer *&p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_V(p,d,pgc);
    HLL(p,d,d->VHs,d->VHn,d->VHe,d->VHw);
    
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

void nhflow_HLLYL::aij_W(lexer *&p,fdm_nhf *&d, int ipol)
{
    double Pval,Qval,wvalbed;
    double dfdx_min, dfdx_plus, dfdy_min, dfdy_plus;
    double detadx,detady;
    

    // HLL W flux sum
    pflux->start_W(p,d,pgc);
    HLL(p,d,d->WHs,d->WHn,d->WHe,d->WHw);
    
    pgc->start1V(p,d->Fx,12);
    pgc->start2V(p,d->Fy,12);
    pgc->start3V(p,d->Fz,12);
    
    LOOP
    WETDRY
    {
    d->FSW[IJK] = ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                +  (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                +  (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
    
    
    // WCALC
    LOOP
    WETDRY
    if(k<p->knoz-1)
    {
    Pval = d->U[IJK];
    Qval = d->V[IJK];
    
    wvalbed=0.0;    
    
        if(k==1)
        {
        dfdx_plus = (d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP];
        dfdx_min  = (d->depth(i,j)-d->depth(i-1,j))/p->DXP[IM1];
    
        detadx = limiter(dfdx_plus,dfdx_min);
        
        dfdy_plus = (d->depth(i,j+1)-d->depth(i,j))/p->DYP[JP];
        dfdy_min  = (d->depth(i,j)-d->depth(i,j-1))/p->DYP[JM1];
    
        detady = limiter(dfdy_plus,dfdy_min);
        
        
        wvalbed =   - Pval*detadx

                    - Qval*detady;
        }
               
        d->test[IJKp1] = d->test[IJK] +  p->WL[IJ]*wvalbed
                            
                    - p->DZN[KP]*(
                            
                    + (d->FEx[IJK] - d->FEx[Im1JK])/p->DXN[IP]  + (d->FEy[IJK] - d->FEy[IJm1K])/p->DYN[JP]*p->y_dir)
                            
                    + p->WL[IJ]*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP])
                            
                    + p->WL[IJ]*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]);
    }

    pgc->start5V(p,d->test,1);
}

void nhflow_HLLYL::aij_E(lexer *&p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_E(p,d,pgc);
    
    HLL_E(p,d);  // -----
    
    LOOP
    WETDRY
    {
    if(p->wet[Ip1J]==0)
    d->FEx[IJK] = 0.0;
    
    if(p->wet[Im1J]==0)
    d->FEx[Im1JK] = 0.0;
    
    if(p->wet[IJp1]==0)
    d->FEy[IJK] = 0.0;
    
    if(p->wet[IJm1]==0)
    d->FEy[IJm1K] = 0.0;
    }
    
    pgc->start1V(p,d->FEx,14);
    pgc->start2V(p,d->FEy,14); 
}

void nhflow_HLLYL::HLL(lexer *&p,fdm_nhf *&d, double *Us, double *Un, double *Ue, double *Uw)
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

void nhflow_HLLYL::HLL_E(lexer *&p, fdm_nhf *&d)
{
    // HLL flux
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        d->FEx[IJK] = d->Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        d->FEx[IJK] = d->Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->FEx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(d->Dn(i,j) - d->Ds(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        d->FEy[IJK] = d->Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        d->FEy[IJK] = d->Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->FEy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(d->Dw(i,j) - d->De(i,j)))/denom;
        }
    }
}


inline double nhflow_HLLYL::limiter(double v1, double v2)
{
    val=0.0;
    
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);

    if(r<0.0)
    phi = 0.0;
    
    if(r>=0.0 && r<0.5)
    phi = 2.0*r;
    
    if(r>=0.5 && r<1.0)
    phi = 1.0;
    
    if(r>=1.0)
    phi = MIN(MIN(r,2.0), 2.0/(1.0+r));
    
    val = 0.5*phi*(v1+v2);

    
    return val;
}

