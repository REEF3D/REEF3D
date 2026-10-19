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
#include"nhflow_solid_flux.h"
#include"nhflow_membrane_beta.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"
#include"nhflow_flux_build_f.h"
#include"vrans.h"
#include"nhflow_flux_face.h"
#include"nhflow_curv.h"

nhflow_HLL::nhflow_HLL (lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_build_f(p,pgc,pBC);
    
    // A 511 2: HLLC (see hllc_sweep_x/y); on the curvilinear grid (G 50) the momentum components are
    // not face-normal/tangential, there the HLL flux is used
    hllc = (p->A511==2 && p->geo_curv!=1) ? 1 : 0;
}

nhflow_HLL::~nhflow_HLL()
{
    delete pflux;
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
// X 15 3 (fb != nullptr): at a face between a cell inside a floating body (FB < 0) and a fluid cell,
// the momentum flux is the physical flux of the fluid side only. The forced body velocity is then
// neither carried nor diffused (HLL jump term) into the water next to the hull.
template<class FL, class FR, class DQ>
inline void hll_sweep_x(lexer *p, fdm_nhf *d, double *F, FL fl, FR fr, DQ dq, const double *fb=nullptr)
{
    int i,j,k;
    
    ULOOP
    {
        if(fb!=nullptr && ((fb[IJK]<0.0) != (fb[Ip1JK]<0.0)))
        {
        F[IJK] = fb[IJK]<0.0 ? fr(i,j,k) : fl(i,j,k);
        continue;
        }
        
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
inline void hll_sweep_y(lexer *p, fdm_nhf *d, double *F, FL fl, FR fr, DQ dq, const double *fb=nullptr)
{
    int i,j,k;
    
    VLOOP
    {
        if(fb!=nullptr && ((fb[IJK]<0.0) != (fb[IJp1K]<0.0)))
        {
        F[IJK] = fb[IJK]<0.0 ? fr(i,j,k) : fl(i,j,k);
        continue;
        }
        
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

// ---------------------------------------------------------------------------
// A 511 2: HLLC flux for the tangential momentum components
//
// For the shallow-water Riemann problem the HLLC flux of D and of the face-normal momentum is
// identical to the HLL flux, provided the contact speed is the consistent one,
//     S* = (SL DR (uR-SR) - SR DL (uL-SL)) / (DR (uR-SR) - DL (uL-SL))   (Toro 2001)
// with D*_L = D*_R = D_HLL. The tangential components phi (V, W across an x-face; U, W across a
// y-face) are carried with the contact:
//     F(phi H) = F_D^HLL * phi_L   if S* >= 0
//                F_D^HLL * phi_R   if S* <  0
// and since F_D^HLL = D* S* with D* > 0, sign(S*) = sign(F_D^HLL): the tangential flux is the HLL
// mass flux of the layer, upwinded with its own sign. It is continuous where the mass flux changes
// sign, and its dissipation is |F_D| dphi instead of SL SR/(SR-SL) dphi ~ c dphi of HLL.
// The normal momentum and the continuity flux stay HLL (hll_sweep_x/y).
//   FL,FR     : physical flux left/right (supersonic branches and floating-body faces, as HLL)
//   PHL,PHR   : tangential velocity left/right, divided by the face porosity squared like FL,FR
// With porosity the mass flux is the HLL flux of the continuity equation n_s deta/dt + div(Q) = 0
// (jump n_s*(D_R - D_L), por_surf_x/y), so the tangential momentum moves with the same volume flux
// as the water: F(phi H) = F_D * phi/n^2 with the face porosity n.
// ---------------------------------------------------------------------------
namespace
{
template<class FL, class FR, class PHL, class PHR>
inline void hllc_sweep_x(lexer *p, fdm_nhf *d, double *F, FL fl, FR fr, PHL phl, PHR phr, const double *fb=nullptr)
{
    int i,j,k;
    
    ULOOP
    {
        if(fb!=nullptr && ((fb[IJK]<0.0) != (fb[Ip1JK]<0.0)))
        {
        F[IJK] = fb[IJK]<0.0 ? fr(i,j,k) : fl(i,j,k);
        continue;
        }
        
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
        
        // HLL mass flux of the layer, the same as the continuity flux (aij_E: jump n_s*(D_R - D_L))
        const double FD = (SR*d->UHs[IJK] - SL*d->UHn[IJK] + SR*SL*nhflow_face::por_surf_x(p,d,i,j)*(d->Dn(i,j) - d->Ds(i,j)))/denom;
        
        F[IJK] = FD>=0.0 ? FD*phl(i,j,k) : FD*phr(i,j,k);
        }
    }
}

template<class FL, class FR, class PHL, class PHR>
inline void hllc_sweep_y(lexer *p, fdm_nhf *d, double *F, FL fl, FR fr, PHL phl, PHR phr, const double *fb=nullptr)
{
    int i,j,k;
    
    VLOOP
    {
        if(fb!=nullptr && ((fb[IJK]<0.0) != (fb[IJp1K]<0.0)))
        {
        F[IJK] = fb[IJK]<0.0 ? fr(i,j,k) : fl(i,j,k);
        continue;
        }
        
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
        
        // HLL mass flux of the layer, the same as the continuity flux (aij_E: jump n_s*(D_R - D_L))
        const double FD = (SR*d->VHe[IJK] - SL*d->VHw[IJK] + SR*SL*nhflow_face::por_surf_y(p,d,i,j)*(d->Dw(i,j) - d->De(i,j)))/denom;
        
        F[IJK] = FD>=0.0 ? FD*phl(i,j,k) : FD*phr(i,j,k);
        }
    }
}
}

namespace
{
// X 15 3: vertical momentum flux at a face between a body cell and a fluid cell from the fluid side only
inline void body_zflux(lexer *p, fdm_nhf *d, const double *fb, const double *Fb, const double *Ft)
{
    if(fb==nullptr)
    return;
    
    int i,j,k;
    
    WLOOP
    if((fb[IJK]<0.0) != (fb[IJKp1]<0.0))
    d->Fz[IJK] = d->omegaF[FIJKp1]*(fb[IJK]<0.0 ? Ft[IJK] : Fb[IJK])/(PORVALNH3*PORVALNH3);
}
}

// curvilinear grid (G 50, nhflow_curv.h): Us/Un, Ve/Vw hold the face-normal velocities; the pressure
// term acts along the face normal, the continuity flux is UH.n
namespace curvface
{
inline double Px(lexer *p, fdm_nhf *d, double eta, int i, int j) { const double g=fabs(p->W22); return 0.5*g*eta*eta + g*eta*d->dfx(i,j); }
inline double Py(lexer *p, fdm_nhf *d, double eta, int i, int j) { const double g=fabs(p->W22); return 0.5*g*eta*eta + g*eta*d->dfy(i,j); }
inline double n1x(lexer *p, int i, int j) { double nx,ny,l; nhflow_curv::normal(p,1,i,j,nx,ny,l); return nx; }
inline double n1y(lexer *p, int i, int j) { double nx,ny,l; nhflow_curv::normal(p,1,i,j,nx,ny,l); return ny; }
inline double n2x(lexer *p, int i, int j) { double nx,ny,l; nhflow_curv::normal(p,2,i,j,nx,ny,l); return nx; }
inline double n2y(lexer *p, int i, int j) { double nx,ny,l; nhflow_curv::normal(p,2,i,j,nx,ny,l); return ny; }
}

void nhflow_HLL::aij_U(lexer *&p,fdm_nhf *&d, int ipol)
{
    // X 15 3: momentum fluxes at body faces from the fluid side only (hll_sweep_x/y, body_zflux)
    const double *fbm = (p->X15==3 && p->X10>0) ? d->FB : nullptr;

    // HLL flux 
    // fused flux build + HLL
    // curvilinear grid (G 50)
    if(p->geo_curv==1)
    {
    hll_sweep_x(p, d, d->Fx,
                [p,d](int i, int j, int k){return d->UHs[IJK]*d->Us[IJK]/(PORVALNH1*PORVALNH1) + curvface::Px(p,d,d->ETAs(i,j),i,j)*curvface::n1x(p,i,j);},
                [p,d](int i, int j, int k){return d->UHn[IJK]*d->Un[IJK]/(PORVALNH1*PORVALNH1) + curvface::Px(p,d,d->ETAn(i,j),i,j)*curvface::n1x(p,i,j);},
                [p,d](int i, int j, int k){return d->UHn[IJK] - d->UHs[IJK];}, fbm);
    
    hll_sweep_y(p, d, d->Fy,
                [p,d](int i, int j, int k){return d->UHe[IJK]*d->Ve[IJK]/(PORVALNH2*PORVALNH2) + curvface::Py(p,d,d->ETAe(i,j),i,j)*curvface::n2x(p,i,j);},
                [p,d](int i, int j, int k){return d->UHw[IJK]*d->Vw[IJK]/(PORVALNH2*PORVALNH2) + curvface::Py(p,d,d->ETAw(i,j),i,j)*curvface::n2x(p,i,j);},
                [p,d](int i, int j, int k){return d->UHw[IJK] - d->UHe[IJK];}, fbm);
    }
    else
    {
        hll_sweep_x(p, d, d->Fx,
                    [p,d](int i, int j, int k){return nhflow_face::U_s(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::U_n(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->UHn[IJK] - d->UHs[IJK];}, fbm);
        
        if(p->j_dir==1 && hllc==1)
        hllc_sweep_y(p, d, d->Fy,
                    [p,d](int i, int j, int k){return nhflow_face::U_e(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::U_w(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->Ue[IJK]/(PORVALNH2*PORVALNH2);},
                    [p,d](int i, int j, int k){return d->Uw[IJK]/(PORVALNH2*PORVALNH2);}, fbm);
        
        else
        if(p->j_dir==1)
        hll_sweep_y(p, d, d->Fy,
                    [p,d](int i, int j, int k){return nhflow_face::U_e(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::U_w(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->UHw[IJK] - d->UHe[IJK];}, fbm);
    }
    
    nhflow_face::zflux(p,d,d->Ub,d->Ut);
    body_zflux(p,d,fbm,d->Ub,d->Ut);
    
    pgc->start1V(p,d->Fx,10);
    pgc->start2V(p,d->Fy,10);
    
    if(p->geo_curv==1)
    nhflow_curv::bc_fluxes(p,d,1,d->Fx,d->Fy);
    pgc->start3V(p,d->Fz,10);
    
    if(phook!=nullptr)
    phook->flux_hook(p,d,hook_id,1,d->Fx,d->Fy);
    
    LOOP
    WETDRY
    {
    if(p->geo_curv==1)
    d->F[IJK] -= (nhflow_curv::divh(p,d->Fx,d->Fy,i,j,k)
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    
    else
    d->F[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

void nhflow_HLL::aij_V(lexer *&p, fdm_nhf *&d, int ipol)
{
    // X 15 3: momentum fluxes at body faces from the fluid side only (hll_sweep_x/y, body_zflux)
    const double *fbm = (p->X15==3 && p->X10>0) ? d->FB : nullptr;

    // HLL flux 
    // fused flux build + HLL
    // curvilinear grid (G 50)
    if(p->geo_curv==1)
    {
    hll_sweep_x(p, d, d->Fx,
                [p,d](int i, int j, int k){return d->VHs[IJK]*d->Us[IJK]/(PORVALNH1*PORVALNH1) + curvface::Px(p,d,d->ETAs(i,j),i,j)*curvface::n1y(p,i,j);},
                [p,d](int i, int j, int k){return d->VHn[IJK]*d->Un[IJK]/(PORVALNH1*PORVALNH1) + curvface::Px(p,d,d->ETAn(i,j),i,j)*curvface::n1y(p,i,j);},
                [p,d](int i, int j, int k){return d->VHn[IJK] - d->VHs[IJK];}, fbm);
    
    hll_sweep_y(p, d, d->Fy,
                [p,d](int i, int j, int k){return d->VHe[IJK]*d->Ve[IJK]/(PORVALNH2*PORVALNH2) + curvface::Py(p,d,d->ETAe(i,j),i,j)*curvface::n2y(p,i,j);},
                [p,d](int i, int j, int k){return d->VHw[IJK]*d->Vw[IJK]/(PORVALNH2*PORVALNH2) + curvface::Py(p,d,d->ETAw(i,j),i,j)*curvface::n2y(p,i,j);},
                [p,d](int i, int j, int k){return d->VHw[IJK] - d->VHe[IJK];}, fbm);
    }
    else
    {
        if(hllc==1)
        hllc_sweep_x(p, d, d->Fx,
                    [p,d](int i, int j, int k){return nhflow_face::V_s(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::V_n(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->Vs[IJK]/(PORVALNH1*PORVALNH1);},
                    [p,d](int i, int j, int k){return d->Vn[IJK]/(PORVALNH1*PORVALNH1);}, fbm);
        
        else
        hll_sweep_x(p, d, d->Fx,
                    [p,d](int i, int j, int k){return nhflow_face::V_s(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::V_n(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->VHn[IJK] - d->VHs[IJK];}, fbm);
        
        if(p->j_dir==1)
        hll_sweep_y(p, d, d->Fy,
                    [p,d](int i, int j, int k){return nhflow_face::V_e(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::V_w(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->VHw[IJK] - d->VHe[IJK];}, fbm);
    }
    
    nhflow_face::zflux(p,d,d->Vb,d->Vt);
    body_zflux(p,d,fbm,d->Vb,d->Vt);
    
    pgc->start1V(p,d->Fx,11);
    pgc->start2V(p,d->Fy,11);
    
    if(p->geo_curv==1)
    nhflow_curv::bc_fluxes(p,d,2,d->Fx,d->Fy);
    pgc->start3V(p,d->Fz,11);
    
    if(phook!=nullptr)
    phook->flux_hook(p,d,hook_id,2,d->Fx,d->Fy);
    
    LOOP
    WETDRY
    {
    if(p->geo_curv==1)
    d->G[IJK] -= (nhflow_curv::divh(p,d->Fx,d->Fy,i,j,k)
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    
    else
    d->G[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

void nhflow_HLL::aij_W(lexer *&p,fdm_nhf *&d, int ipol)
{
    // X 15 3: momentum fluxes at body faces from the fluid side only (hll_sweep_x/y, body_zflux)
    const double *fbm = (p->X15==3 && p->X10>0) ? d->FB : nullptr;

    // HLL flux 
    // fused flux build + HLL
    // curvilinear grid (G 50)
    if(p->geo_curv==1)
    {
    hll_sweep_x(p, d, d->Fx,
                [p,d](int i, int j, int k){return d->WHs[IJK]*d->Us[IJK]/(PORVALNH1*PORVALNH1);},
                [p,d](int i, int j, int k){return d->WHn[IJK]*d->Un[IJK]/(PORVALNH1*PORVALNH1);},
                [p,d](int i, int j, int k){return d->WHn[IJK] - d->WHs[IJK];}, fbm);
    
    hll_sweep_y(p, d, d->Fy,
                [p,d](int i, int j, int k){return d->WHe[IJK]*d->Ve[IJK]/(PORVALNH2*PORVALNH2);},
                [p,d](int i, int j, int k){return d->WHw[IJK]*d->Vw[IJK]/(PORVALNH2*PORVALNH2);},
                [p,d](int i, int j, int k){return d->WHw[IJK] - d->WHe[IJK];}, fbm);
    }
    else
    {
        if(hllc==1)
        {
        hllc_sweep_x(p, d, d->Fx,
                    [p,d](int i, int j, int k){return nhflow_face::W_s(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::W_n(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->Ws[IJK]/(PORVALNH1*PORVALNH1);},
                    [p,d](int i, int j, int k){return d->Wn[IJK]/(PORVALNH1*PORVALNH1);}, fbm);
        
        if(p->j_dir==1)
        hllc_sweep_y(p, d, d->Fy,
                    [p,d](int i, int j, int k){return nhflow_face::W_e(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::W_w(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->We[IJK]/(PORVALNH2*PORVALNH2);},
                    [p,d](int i, int j, int k){return d->Ww[IJK]/(PORVALNH2*PORVALNH2);}, fbm);
        }
        
        else
        {
        hll_sweep_x(p, d, d->Fx,
                    [p,d](int i, int j, int k){return nhflow_face::W_s(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::W_n(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->WHn[IJK] - d->WHs[IJK];}, fbm);
        
        if(p->j_dir==1)
        hll_sweep_y(p, d, d->Fy,
                    [p,d](int i, int j, int k){return nhflow_face::W_e(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return nhflow_face::W_w(p,d,i,j,k);},
                    [p,d](int i, int j, int k){return d->WHw[IJK] - d->WHe[IJK];}, fbm);
        }
    }
    
    nhflow_face::zflux(p,d,d->Wb,d->Wt);
    body_zflux(p,d,fbm,d->Wb,d->Wt);
    
    pgc->start1V(p,d->Fx,12);
    pgc->start2V(p,d->Fy,12);
    
    if(p->geo_curv==1)
    nhflow_curv::bc_fluxes(p,d,3,d->Fx,d->Fy);
    pgc->start3V(p,d->Fz,12);
    
    if(phook!=nullptr)
    phook->flux_hook(p,d,hook_id,3,d->Fx,d->Fy);
    
    LOOP
    WETDRY
    {
    if(p->geo_curv==1)
    d->H[IJK] -= (nhflow_curv::divh(p,d->Fx,d->Fy,i,j,k)
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    
    else
    d->H[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

void nhflow_HLL::aij_E(lexer *&p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    // fused flux build + HLL (continuity: F = UH, VH;  q = D)
    // curvilinear grid (G 50): the continuity flux UH.n
    if(p->geo_curv==1)
    {
    hll_sweep_x(p, d, d->FEx,
                [p,d](int i, int j, int k){return d->UHs[IJK]*curvface::n1x(p,i,j) + d->VHs[IJK]*curvface::n1y(p,i,j);},
                [p,d](int i, int j, int k){return d->UHn[IJK]*curvface::n1x(p,i,j) + d->VHn[IJK]*curvface::n1y(p,i,j);},
                [p,d](int i, int j, int k){return nhflow_face::por_surf_x(p,d,i,j)*(d->Dn(i,j) - d->Ds(i,j));});
    
    hll_sweep_y(p, d, d->FEy,
                [p,d](int i, int j, int k){return d->UHe[IJK]*curvface::n2x(p,i,j) + d->VHe[IJK]*curvface::n2y(p,i,j);},
                [p,d](int i, int j, int k){return d->UHw[IJK]*curvface::n2x(p,i,j) + d->VHw[IJK]*curvface::n2y(p,i,j);},
                [p,d](int i, int j, int k){return nhflow_face::por_surf_y(p,d,i,j)*(d->Dw(i,j) - d->De(i,j));});
    }
    else
    {
    hll_sweep_x(p, d, d->FEx,
                [p,d](int i, int j, int k){return d->UHs[IJK];},
                [p,d](int i, int j, int k){return d->UHn[IJK];},
                [p,d](int i, int j, int k){return nhflow_face::por_surf_x(p,d,i,j)*(d->Dn(i,j) - d->Ds(i,j));});
    
    if(p->j_dir==1)
    hll_sweep_y(p, d, d->FEy,
                [p,d](int i, int j, int k){return d->VHe[IJK];},
                [p,d](int i, int j, int k){return d->VHw[IJK];},
                [p,d](int i, int j, int k){return nhflow_face::por_surf_y(p,d,i,j)*(d->Dw(i,j) - d->De(i,j));});
    }
    
    // membranes (X 330): Rhie-Chow continuity flux at faces next to membrane cells, see nhflow_membrane_beta.h
    if(d->MBETA!=nullptr)
    nhflow_membrane_rc_flux(p,d);
    
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
    
    // cells inside FEM structures (Z 30) and direct forcing solids: no continuity flux through their faces
    if(d->solid_flux>0)
    nhflow_solid_flux(p,d);
    
    pgc->start1V(p,d->FEx,14);
    pgc->start2V(p,d->FEy,14); 
    
    if(p->geo_curv==1)
    nhflow_curv::bc_fluxes(p,d,4,d->FEx,d->FEy);
    
    if(phook!=nullptr)
    phook->flux_hook(p,d,hook_id,4,d->FEx,d->FEy);
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
        
        d->FEx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*nhflow_face::por_surf_x(p,d,i,j)*(d->Dn(i,j) - d->Ds(i,j)))/denom;
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
        
        d->FEy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*nhflow_face::por_surf_y(p,d,i,j)*(d->Dw(i,j) - d->De(i,j)))/denom;
        }
    }
}
