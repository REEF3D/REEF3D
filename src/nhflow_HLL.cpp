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

#include"nhflow_HLL.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"
#include"nhflow_flux_build_f.h"
#include"vrans.h"
#include"nhflow_flux_face.h"

nhflow_HLL::nhflow_HLL (lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_build_f(p,pgc,pBC);
}

nhflow_HLL::~nhflow_HLL()
{
}

void nhflow_HLL::precalc(lexer* p, fdm_nhf* d, int ipolL, slice &eta)
{
}

void nhflow_HLL::start(lexer *&p, fdm_nhf *&d, int ipol, slice &eta, double *FH)
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

// ---------------------------------------------------------------------------
// Fused flux build + HLL: the left/right physical fluxes of a face are
// evaluated on the fly (nhflow_flux_face.h, the same expressions that
// nhflow_flux_build_f stores in d->Fs/Fn/Fe/Fw) instead of being written to
// and re-read from four full 3D arrays. Only the fluxes of the branch that is
// taken are evaluated. Loop ranges and arithmetic are those of
// flux_build_f::start_* followed by HLL()/HLL_E().
//   FL,FR : physical flux left/right of the face
//   DQ    : (qR - qL), the jump of the conserved variable
// ---------------------------------------------------------------------------
namespace
{
template<class FL, class FR, class DQ>
inline void hll_sweep_x(lexer *p, fdm_nhf *d, double *F, FL fl, FR fr, DQ dq)
{
    int i,j,k;
    
    ULOOP
    {
        const double SL = d->Ss[IJK];
        const double SR = d->Sn[IJK];
        
        if(SL>=0.0)
        F[IJK] = fl(i,j,k);
        
        else
        if(SR<=0.0)
        F[IJK] = fr(i,j,k);
        
        else
        {
        double denom = SR-SL;
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        F[IJK] = (SR*fl(i,j,k) - SL*fr(i,j,k) + SR*SL*dq(i,j,k))/denom;
        }
    }
}

template<class FL, class FR, class DQ>
inline void hll_sweep_y(lexer *p, fdm_nhf *d, double *F, FL fl, FR fr, DQ dq)
{
    int i,j,k;
    
    VLOOP
    {
        const double SL = d->Se[IJK];
        const double SR = d->Sw[IJK];
        
        if(SL>=0.0)
        F[IJK] = fl(i,j,k);
        
        else
        if(SR<=0.0)
        F[IJK] = fr(i,j,k);
        
        else
        {
        double denom = SR-SL;
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        F[IJK] = (SR*fl(i,j,k) - SL*fr(i,j,k) + SR*SL*dq(i,j,k))/denom;
        }
    }
}
}

void nhflow_HLL::aij_U(lexer *&p,fdm_nhf *&d, int ipol)
{
    // HLL flux 
    // fused flux build + HLL
    hll_sweep_x(p, d, d->Fx,
                [p,d](int i, int j, int k){return nhflow_face::U_s(p,d,i,j,k);},
                [p,d](int i, int j, int k){return nhflow_face::U_n(p,d,i,j,k);},
                [p,d](int i, int j, int k){return d->UHn[IJK] - d->UHs[IJK];});
    
    if(p->j_dir==1)
    hll_sweep_y(p, d, d->Fy,
                [p,d](int i, int j, int k){return nhflow_face::U_e(p,d,i,j,k);},
                [p,d](int i, int j, int k){return nhflow_face::U_w(p,d,i,j,k);},
                [p,d](int i, int j, int k){return d->UHw[IJK] - d->UHe[IJK];});
    
    nhflow_face::zflux(p,d,d->Ub,d->Ut);
    
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

void nhflow_HLL::aij_V(lexer *&p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    // fused flux build + HLL
    hll_sweep_x(p, d, d->Fx,
                [p,d](int i, int j, int k){return nhflow_face::V_s(p,d,i,j,k);},
                [p,d](int i, int j, int k){return nhflow_face::V_n(p,d,i,j,k);},
                [p,d](int i, int j, int k){return d->VHn[IJK] - d->VHs[IJK];});
    
    if(p->j_dir==1)
    hll_sweep_y(p, d, d->Fy,
                [p,d](int i, int j, int k){return nhflow_face::V_e(p,d,i,j,k);},
                [p,d](int i, int j, int k){return nhflow_face::V_w(p,d,i,j,k);},
                [p,d](int i, int j, int k){return d->VHw[IJK] - d->VHe[IJK];});
    
    nhflow_face::zflux(p,d,d->Vb,d->Vt);
    
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

void nhflow_HLL::aij_W(lexer *&p,fdm_nhf *&d, int ipol)
{
    // HLL flux 
    // fused flux build + HLL
    hll_sweep_x(p, d, d->Fx,
                [p,d](int i, int j, int k){return nhflow_face::W_s(p,d,i,j,k);},
                [p,d](int i, int j, int k){return nhflow_face::W_n(p,d,i,j,k);},
                [p,d](int i, int j, int k){return d->WHn[IJK] - d->WHs[IJK];});
    
    if(p->j_dir==1)
    hll_sweep_y(p, d, d->Fy,
                [p,d](int i, int j, int k){return nhflow_face::W_e(p,d,i,j,k);},
                [p,d](int i, int j, int k){return nhflow_face::W_w(p,d,i,j,k);},
                [p,d](int i, int j, int k){return d->WHw[IJK] - d->WHe[IJK];});
    
    nhflow_face::zflux(p,d,d->Wb,d->Wt);
    
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

void nhflow_HLL::aij_E(lexer *&p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    // fused flux build + HLL (continuity: F = UH, VH;  q = D)
    hll_sweep_x(p, d, d->FEx,
                [p,d](int i, int j, int k){return d->UHs[IJK];},
                [p,d](int i, int j, int k){return d->UHn[IJK];},
                [p,d](int i, int j, int k){return d->Dn(i,j) - d->Ds(i,j);});
    
    if(p->j_dir==1)
    hll_sweep_y(p, d, d->FEy,
                [p,d](int i, int j, int k){return d->VHe[IJK];},
                [p,d](int i, int j, int k){return d->VHw[IJK];},
                [p,d](int i, int j, int k){return d->Dw(i,j) - d->De(i,j);});
    
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

void nhflow_HLL::HLL(lexer *&p,fdm_nhf *&d, double *Us, double *Un, double *Ue, double *Uw)
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

void nhflow_HLL::HLL_E(lexer *&p, fdm_nhf *&d)
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
