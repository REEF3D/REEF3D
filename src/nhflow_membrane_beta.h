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


#ifndef NHFLOW_MEMBRANE_BETA_H_
#define NHFLOW_MEMBRANE_BETA_H_

// Membrane mobility (X 330) in the NHFLOW projection.
//
// The membrane forcing is integrated implicitly, (1 + a K H)(u - u_m) = (u* - u_m) - a/rho grad(q),
// so in the smeared membrane layer the pressure acts with the reduced mobility beta = 1/(1 + a K H)
// (d->MBETA, cell centred). The mobility enters the projection in three consistent places:
//
// 1. Poisson matrix (nhflow_membrane_row). Row n belongs to the pressure node (i,j,k) at the bottom
//    face of cell k, between cell k-1 (below) and cell k (above):
//      horizontal faces   harmonic mean of the two cells' beta, averaged over cells k-1, k with the
//                         weights of the velocity average in the divergence
//      vertical faces     beta of the cell between the nodes
//      sigma cross terms  node average bF (so that a constant pressure stays in the null space)
//
// 2. Velocity correction (nhflow_membrane_gradx/grady). The collocated correction uses the wide
//    gradient (P_i+1 - P_i-1)/(2 dx). Written as the distance-weighted mean of the two compact face
//    gradients, the mobility is applied per face:
//      dU_i = -a/rho [ w+ beta_i+1/2 (P_i+1 - P_i)/dx+ + w- beta_i-1/2 (P_i - P_i-1)/dx- ]
//    Without membrane (beta = 1) this is identical to the wide gradient. With the cell value instead,
//    a cell next to the layer is corrected with mobility 1 by the pressure difference across the
//    layer that the matrix only sees through the small face mobility: the projection over-corrects
//    and diverges in 3D and with the incremental scheme (A 520 2).
//
// 3. Continuity flux (nhflow_membrane_rc_flux). Rhie-Chow face velocity near the membrane,
//      U_f = (U_i + U_i+1)/2 + (dU_i + dU_i+1)/2 - dU_f,   dU_f = -a/rho beta_f (P_i+1 - P_i)/dx
//    i.e. the face velocity the compact Poisson matrix makes divergence free. The free surface then
//    follows the projected face field: no mass moves through the membrane layer or below the bag floor
//    that the projection does not see. dU_f of the last projection is kept in d->MRCX, d->MRCY.

#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"increment.h"
#include"vrans_definitions.h"
#include"nhflow_thinbody.h"

// mobility of the vertical velocity correction in cell (i,j,k): compact difference of the nodes k, k+1,
// exactly the vertical face of the Poisson matrix
#define MBETAVAL (d->MBZ!=nullptr ? d->MBZ[IJK] : (d->MBETA!=nullptr ? d->MBETA[IJK] : 1.0))


inline double nhflow_membrane_harm(double a, double b)
{
    return 2.0*a*b/(a+b);
}

// Link mobilities. Layer mode (default): harmonic mean of the cell mobilities beta of the two cells. Link mode
// (membrane.dat 'mobility link', d->MBX/MBY/MBZ allocated): the mobility of the link itself, small only where the
// link between the two cell centres (nodes for the vertical links) crosses the membrane; the fluid in the layer
// moves freely along the membrane (anisotropic: normal flow blocked on the crossing links, tangential flow free).
// q: cell (i,j,k), qp: its +x / +y neighbour
inline double nhflow_mbx(const fdm_nhf *d, int q, int qp)
{
    return d->MBX!=nullptr ? d->MBX[q] : nhflow_membrane_harm(d->MBETA[q],d->MBETA[qp]);
}

inline double nhflow_mby(const fdm_nhf *d, int q, int qp)
{
    return d->MBY!=nullptr ? d->MBY[q] : nhflow_membrane_harm(d->MBETA[q],d->MBETA[qp]);
}

inline double nhflow_mbz(const fdm_nhf *d, int q)
{
    return d->MBZ!=nullptr ? d->MBZ[q] : d->MBETA[q];
}

// cell whose horizontal correction differs from the plain wide gradient
inline bool nhflow_membrane_active(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    const double *B = d->MBETA;
    
    if(B[IJK]<1.0 || B[Ip1JK]<1.0 || B[Im1JK]<1.0)
    return true;
    
    if(p->j_dir==1 && (B[IJp1K]<1.0 || B[IJm1K]<1.0))
    return true;
    
    return false;
}

// horizontal pressure gradient of the velocity correction with face mobilities (without the sigma term);
// P at the nodes, averaged to the cell level of layer k as in nhflow_pjm::ucorr
inline double nhflow_membrane_gradx(lexer *p, fdm_nhf *d, const double *P, int i, int j, int k)
{
    const int marge = increment::marge;
    const double *B = d->MBETA;
    
    const nhflow_thinbody *tb = d->thinbody;
    const double Pc = tb!=nullptr ? tb->pcell(p,P,i,j,k)   : 0.5*(P[FIJK]+P[FIJKp1]);
    const double Pn = tb!=nullptr ? tb->pcell(p,P,i+1,j,k) : 0.5*(P[FIp1JK]+P[FIp1JKp1]);
    const double Ps = tb!=nullptr ? tb->pcell(p,P,i-1,j,k) : 0.5*(P[FIm1JK]+P[FIm1JKp1]);
    
    const double dxp = p->DXP[IP], dxm = p->DXP[IM1];
    
    return (nhflow_mbx(d,IJK,Ip1JK)*(Pn-Pc) + nhflow_mbx(d,Im1JK,IJK)*(Pc-Ps))/(dxp+dxm);
}

inline double nhflow_membrane_grady(lexer *p, fdm_nhf *d, const double *P, int i, int j, int k)
{
    const int marge = increment::marge;
    const double *B = d->MBETA;
    
    const nhflow_thinbody *tb = d->thinbody;
    const double Pc = tb!=nullptr ? tb->pcell(p,P,i,j,k)   : 0.5*(P[FIJK]+P[FIJKp1]);
    const double Pw = tb!=nullptr ? tb->pcell(p,P,i,j+1,k) : 0.5*(P[FIJp1K]+P[FIJp1Kp1]);
    const double Pe = tb!=nullptr ? tb->pcell(p,P,i,j-1,k) : 0.5*(P[FIJm1K]+P[FIJm1Kp1]);
    
    const double dyp = p->DYP[JP], dym = p->DYP[JM1];
    
    return (nhflow_mby(d,IJK,IJp1K)*(Pw-Pc) + nhflow_mby(d,IJm1K,IJK)*(Pc-Pe))/(dyp+dym);
}

inline void nhflow_membrane_row(lexer *p, fdm_nhf *d, int i, int j, int k, int n, double ct, double cb, double sxx, double rhs0)
{
    // ct, cb: vertical Laplacian coefficients of the row (M.t, M.b without the sigxx term)
    // sxx:    sigxx first-derivative coefficient (M.t = ct - sxx, M.b = cb + sxx)
    //
    // The horizontal face mobility of the node is the average of the face mobilities of the two cells
    // it sits between, weighted like the vertical average of the velocities in the divergence
    // (nhflow_pjm::rhs). The matrix then matches divergence o correction (nhflow_membrane_gradx) for
    // smooth pressure fields; a min over the two cells makes the matrix weaker than the correction
    // where beta changes vertically (bag floor), and the incremental scheme diverges there.
    const int marge = increment::marge;
    const double *B = d->MBETA;
    
    const double fac = k>0 ? p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]) : 0.0;
    
    double bc, bcm, bn, bs, bw, be;
    
    if(d->MBX!=nullptr)
    {
        // link mode: the node links are the cell links of cells k and k-1, weighted like the velocities
        bc  = d->MBZ[IJK];
        bcm = k>0 ? d->MBZ[IJKm1] : bc;
        
        bn = (1.0-fac)*d->MBX[IJK]   + fac*(k>0 ? d->MBX[IJKm1]   : d->MBX[IJK]);
        bs = (1.0-fac)*d->MBX[Im1JK] + fac*(k>0 ? d->MBX[Im1JKm1] : d->MBX[Im1JK]);
        bw = p->j_dir==1 ? (1.0-fac)*d->MBY[IJK]   + fac*(k>0 ? d->MBY[IJKm1]   : d->MBY[IJK])   : 1.0;
        be = p->j_dir==1 ? (1.0-fac)*d->MBY[IJm1K] + fac*(k>0 ? d->MBY[IJm1Km1] : d->MBY[IJm1K]) : 1.0;
        
        // sharp thin bodies: a horizontal node link whose nodes lie on opposite sides of the body is blocked as well
        // (next to a sloping floor the node control volume straddles the floor, the cell layers do not see it), and the
        // half face of a layer whose two cells lie on the other side of the body than the node is closed
        if(d->thinbody!=nullptr)
        {
            const nhflow_thinbody *tb = d->thinbody;
            const double bmin = 1.0e-4;
            const double sN = tb->nside(p,i,j,k);
            
            auto lay = [&](const double *MB, int qa, int qb) {return tb->other_side(qa,qb,sN) ? MIN(MB[qa],bmin) : MB[qa];};
            
            bn = (1.0-fac)*lay(d->MBX,IJK,Ip1JK)   + fac*(k>0 ? lay(d->MBX,IJKm1,Ip1JKm1)   : lay(d->MBX,IJK,Ip1JK));
            bs = (1.0-fac)*lay(d->MBX,Im1JK,IJK)   + fac*(k>0 ? lay(d->MBX,Im1JKm1,IJKm1)   : lay(d->MBX,Im1JK,IJK));
            
            if(p->j_dir==1)
            {
            bw = (1.0-fac)*lay(d->MBY,IJK,IJp1K)   + fac*(k>0 ? lay(d->MBY,IJKm1,IJp1Km1)   : lay(d->MBY,IJK,IJp1K));
            be = (1.0-fac)*lay(d->MBY,IJm1K,IJK)   + fac*(k>0 ? lay(d->MBY,IJm1Km1,IJKm1)   : lay(d->MBY,IJm1K,IJK));
            }
            
            if(tb->node_cut_x(p,i,j,k))   bn = MIN(bn, MIN(MIN(d->MBX[IJK],k>0?d->MBX[IJKm1]:1.0),bmin));
            if(tb->node_cut_x(p,i-1,j,k)) bs = MIN(bs, MIN(MIN(d->MBX[Im1JK],k>0?d->MBX[Im1JKm1]:1.0),bmin));
            
            if(p->j_dir==1)
            {
            if(tb->node_cut_y(p,i,j,k))   bw = MIN(bw, MIN(MIN(d->MBY[IJK],k>0?d->MBY[IJKm1]:1.0),bmin));
            if(tb->node_cut_y(p,i,j-1,k)) be = MIN(be, MIN(MIN(d->MBY[IJm1K],k>0?d->MBY[IJm1Km1]:1.0),bmin));
            }
        }
    }
    else
    {
        bc  = B[IJK];
        bcm = k>0 ? B[IJKm1] : bc;
        
        auto hmean = [&](double a0, double a1, double b0, double b1)
        {
            return (1.0-fac)*nhflow_membrane_harm(a0,a1) + fac*nhflow_membrane_harm(b0,b1);
        };
        
        bn = hmean(bc,B[Ip1JK], bcm,k>0?B[Ip1JKm1]:B[Ip1JK]);
        bs = hmean(bc,B[Im1JK], bcm,k>0?B[Im1JKm1]:B[Im1JK]);
        bw = p->j_dir==1 ? hmean(bc,B[IJp1K], bcm,k>0?B[IJp1Km1]:B[IJp1K]) : 1.0;
        be = p->j_dir==1 ? hmean(bc,B[IJm1K], bcm,k>0?B[IJm1Km1]:B[IJm1K]) : 1.0;
    }
    const double bF = (1.0-fac)*bc + fac*bcm;
    const double bt = bc;
    const double bb = bcm;
    
    if(bn==1.0 && bs==1.0 && bw==1.0 && be==1.0 && bt==1.0 && bb==1.0)
    return;
    
    d->M.n[n] *= bn;
    d->M.s[n] *= bs;
    d->M.w[n] *= bw;
    d->M.e[n] *= be;
    
    d->M.t[n] = bt*ct - bF*sxx;
    d->M.b[n] = bb*cb + bF*sxx;
    
    d->M.p[n] = -(d->M.n[n] + d->M.s[n] + d->M.w[n] + d->M.e[n]) - bt*ct - bb*cb;
    
    // sharp thin bodies: the explicit sigma cross terms of a node next to a blocked link would reach across the body
    const double bX = d->thinbody!=nullptr ? MIN(bF,MIN(MIN(MIN(bn,bs),MIN(bw,be)),MIN(bt,bb))) : bF;
    
    d->rhsvec.V[n] = rhs0 + bX*(d->rhsvec.V[n] - rhs0);
    
    // sharp thin bodies: a node control volume enclosed by the body on all sides (a fold of a flexible membrane, a
    // corner pocket) has no open link; its divergence is scaled with the largest link mobility, so that the near
    // singular row does not turn the wall fluxes of the pocket into a pressure spike (no change for open nodes)
    if(d->thinbody!=nullptr)
    {
        const double bmx = MAX(MAX(MAX(bn,bs),MAX(bw,be)),MAX(bt,bb));
        
        if(bmx<1.0)
        d->rhsvec.V[n] *= bmx;
    }
}

// face correction velocities dU_f of the projection with the total pressure P (after the correction),
// a = alpha dt, for the Rhie-Chow continuity flux of the next stage
// VRANS: scaled by the face value of CPORNH, as the Poisson matrix and ucorr/vcorr (= 1 without porosity)
inline void nhflow_membrane_rc_store(lexer *p, fdm_nhf *d, ghostcell *pgc, const double *P, double a)
{
    int i,j,k;
    const int marge = increment::marge;
    const double *B = d->MBETA;
    
    LOOP
    {
        const nhflow_thinbody *tb = d->thinbody;
        const double Pc = tb!=nullptr ? tb->pcell(p,P,i,j,k)   : 0.5*(P[FIJK]+P[FIJKp1]);
        const double Pn = tb!=nullptr ? tb->pcell(p,P,i+1,j,k) : 0.5*(P[FIp1JK]+P[FIp1JKp1]);
        
        const double cx = 0.5*(CPORNHval(d->POR[IJK]) + CPORNHval(d->POR[Ip1JK]));
        
        d->MRCX[IJK] = -a*cx/p->W1*nhflow_mbx(d,IJK,Ip1JK)*(Pn-Pc)/p->DXP[IP];
        
        if(p->j_dir==1)
        {
        const double Pw = tb!=nullptr ? tb->pcell(p,P,i,j+1,k) : 0.5*(P[FIJp1K]+P[FIJp1Kp1]);
        const double cy = 0.5*(CPORNHval(d->POR[IJK]) + CPORNHval(d->POR[IJp1K]));
        
        d->MRCY[IJK] = -a*cy/p->W1*nhflow_mby(d,IJK,IJp1K)*(Pw-Pc)/p->DYP[JP];
        }
    }
    
    pgc->start4V(p,d->MRCX,1);
    
    if(p->j_dir==1)
    pgc->start4V(p,d->MRCY,1);
}

// U = UH/WL after a correction, for a further projection pass within the stage
inline void nhflow_membrane_velupdate(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *UH, double *VH, double *WH)
{
    int i,j,k;
    
    LOOP
    if(p->wet[IJ]==1)
    {
        const double wl = fabs(WL(i,j))>1.0e-20 ? WL(i,j) : 1.0e20;
        d->U[IJK] = UH[IJK]/wl;
        d->V[IJK] = VH[IJK]/wl;
        d->W[IJK] = WH[IJK]/wl;
    }
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
}

inline bool nhflow_membrane_face_cell(fdm_nhf *d, int q)
{
    return d->MBETA[q]<1.0 || d->MCHI[q]>0.0;
}

// Rhie-Chow continuity flux at faces within two cells of membrane cells (beta < 1) or of the region
// below the bag floor with the prescribed static pressure (chi > 0)
inline void nhflow_membrane_rc_flux(lexer *p, fdm_nhf *d)
{
    int i,j,k;
    const int marge = increment::marge;
    
    ULOOP
    if(nhflow_membrane_face_cell(d,IJK) || nhflow_membrane_face_cell(d,Ip1JK) || nhflow_membrane_face_cell(d,Im1JK) || nhflow_membrane_face_cell(d,Ip2JK))
    {
        const double dxs = p->DXP[IM1], dxc = p->DXP[IP], dxn = p->DXP[IP1];
        
        // cell corrections as distance-weighted means of the face corrections
        const double dUi  = (dxc*d->MRCX[IJK]   + dxs*d->MRCX[Im1JK])/(dxc+dxs);
        const double dUi1 = (dxn*d->MRCX[Ip1JK] + dxc*d->MRCX[IJK])/(dxn+dxc);
        
        // link mode: U_f = avg(U*) + dU_f = avg(U) - avg(dU) + dU_f (the cells next to a blocked link are corrected
        // through their open faces, the blocked face is not)
        const double rc = d->MPROJ>1 ? 0.0 : (d->MBX!=nullptr ? d->MRCX[IJK] - 0.5*(dUi + dUi1) : 0.5*(dUi + dUi1) - d->MRCX[IJK]);
        const double Uf = 0.5*(d->U[IJK] + d->U[Ip1JK]) + rc;
        
        d->FEx[IJK] = 0.5*(d->Ds(i,j) + d->Dn(i,j))*Uf;
    }
    
    if(p->j_dir==1)
    VLOOP
    if(nhflow_membrane_face_cell(d,IJK) || nhflow_membrane_face_cell(d,IJp1K) || nhflow_membrane_face_cell(d,IJm1K) || nhflow_membrane_face_cell(d,IJp2K))
    {
        const double dys = p->DYP[JM1], dyc = p->DYP[JP], dyn = p->DYP[JP1];
        
        const double dVj  = (dyc*d->MRCY[IJK]   + dys*d->MRCY[IJm1K])/(dyc+dys);
        const double dVj1 = (dyn*d->MRCY[IJp1K] + dyc*d->MRCY[IJK])/(dyn+dyc);
        
        const double rc = d->MPROJ>1 ? 0.0 : (d->MBX!=nullptr ? d->MRCY[IJK] - 0.5*(dVj + dVj1) : 0.5*(dVj + dVj1) - d->MRCY[IJK]);
        const double Vf = 0.5*(d->V[IJK] + d->V[IJp1K]) + rc;
        
        d->FEy[IJK] = 0.5*(d->De(i,j) + d->Dw(i,j))*Vf;
    }
}

#endif
