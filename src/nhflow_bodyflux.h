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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_BODYFLUX_H_
#define NHFLOW_BODYFLUX_H_

// NHFLOW direct forcing with the body in the projection (X 15 11): the fluxes at the faces of the body.
//
// A face with the open fraction b (d->MBX, d->MBY; the vertical link through a cell d->MBZ, 1 - the fraction of the
// cell inside the body; sixdof_obj_nhflow::forcing_kernel_fraction) carries the water through its open part and the
// body through the rest:
//
//      u_face = b (U_i + U_i+1)/2 + (1 - b) u_b
//
// with u_b the normal velocity of the body at the face (d->FBUX, d->FBUY, vertically d->FBWZ at the cell centres).
// The same face velocities are used
//   - in the divergence of the pressure projection (rhs, called by nhflow_pjm::rhs and nhflow_pjm_corr::rhs): the
//     projection makes the water flow around the hull with the velocity of the hull at its faces, whatever the cell
//     velocities next to it are, and no pressure is spent on water that the closed faces cannot let through,
//   - in the continuity flux of the water level (nhflow_HLL::aij_E): no water moves through the hull, the columns
//     next to it take the volume it displaces,
//   - and, with the same weights, in the momentum fluxes (nhflow_HLL, hll_sweep_x/y): b F_HLL + (1 - b) F_wall.
//
// Near the body (d->FBNEAR, the columns within nnear cells of a column with a cell of the body) the face velocity
// of the open part gets the Rhie-Chow correction of the last projection (rcx, rcy): the mean of two cell velocities
// does not see a pressure that alternates from column to column, which the faces of the body excite; such a
// pressure then stays and drives the vertical velocity in alternate columns up and down.

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"

namespace nhflow_bodyflux
{
    const int nnear = 4;
    
    inline bool active(lexer *p, const fdm_nhf *d)
    {
        return p->X15==11 && p->X10>0 && d->FBUX!=nullptr && d->MBX!=nullptr;
    }
    
    // Rhie-Chow correction of the velocity at the x face (i+1/2,j,k) / y face (i,j+1/2,k): the correction of the face
    // by the pressure of the last projection (d->MRCX, MRCY, nhflow_membrane_rc_store) minus the mean of the corrections
    // of its two cells (each the mean of its two faces, as the pressure gradient of the cell)
    inline double rcx(lexer *p, const fdm_nhf *d, int i, int j, int k)
    {
        const int marge = increment::marge;
        
        if(d->FBNEAR[IJ]==0 && d->FBNEAR[Ip1J]==0)
        return 0.0;
        
        const double dxs = p->DXP[IM1], dxc = p->DXP[IP], dxn = p->DXP[IP1];
        const double dUi  = (dxc*d->MRCX[IJK]   + dxs*d->MRCX[Im1JK])/(dxc+dxs);
        const double dUi1 = (dxn*d->MRCX[Ip1JK] + dxc*d->MRCX[IJK])/(dxn+dxc);
        
        return d->MRCX[IJK] - 0.5*(dUi + dUi1);
    }
    
    inline double rcy(lexer *p, const fdm_nhf *d, int i, int j, int k)
    {
        const int marge = increment::marge;
        
        if(d->FBNEAR[IJ]==0 && d->FBNEAR[IJp1]==0)
        return 0.0;
        
        const double dys = p->DYP[JM1], dyc = p->DYP[JP], dyn = p->DYP[JP1];
        const double dVj  = (dyc*d->MRCY[IJK]   + dys*d->MRCY[IJm1K])/(dyc+dys);
        const double dVj1 = (dyn*d->MRCY[IJp1K] + dyc*d->MRCY[IJK])/(dyn+dyc);
        
        return d->MRCY[IJK] - 0.5*(dVj + dVj1);
    }
    
    // Right-hand side of the Poisson equation (-div/(a dt), a dt = adt) for the cell velocities U, V, W: the change of
    // the divergence of the mean cell velocities to that of the face velocities. Row n is the pressure node at the
    // bottom of cell k: it averages the x / y faces of the layers k (weight 1 - fac) and k-1 (fac), as the horizontal
    // divergence of nhflow_pjm::rhs, vertically the links through the cells k and k-1.
    // ubx, uby, wbz: the velocity of the body at the faces (d->FBUX, ... or those of a unit motion, added mass)
    // rc: with the Rhie-Chow correction (not for the velocity field of a unit motion)
    inline void rhs(lexer *p, fdm_nhf *d, const double *U, const double *V, const double *W, 
                    const double *ubx, const double *uby, const double *wbz, double adt, bool rc)
    {
        int i,j,k;
        const int marge = increment::marge;
        const int *nb = d->FBNEAR;
        
        auto dx = [&](int qc, int qn, int qf, int di, int dk) 
        {return (1.0 - d->MBX[qf])*(ubx[qf] - 0.5*(U[qc] + U[qn])) + (rc ? d->MBX[qf]*rcx(p,d,i+di,j,k+dk) : 0.0);};
        
        auto dy = [&](int qc, int qn, int qf, int dj, int dk) 
        {return (1.0 - d->MBY[qf])*(uby[qf] - 0.5*(V[qc] + V[qn])) + (rc ? d->MBY[qf]*rcy(p,d,i,j+dj,k+dk) : 0.0);};
        
        auto dz = [&](int qc) {return (1.0 - d->MBZ[qc])*(wbz[qc] - W[qc]);};
        
        int n=0;
        LOOP
        {
            if(p->wet[IJ]==1 && p->deep[IJ]==1 && p->flag7[FIJK]>0 
               && (nb[IJ]==1 || nb[Ip1J]==1 || nb[Im1J]==1 || nb[IJp1]==1 || nb[IJm1]==1))
            {
                const bool bed = (k==0);
                const double fac = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
                
                double div = (1.0-fac)*(dx(IJK,Ip1JK,IJK,0,0) - dx(Im1JK,IJK,Im1JK,-1,0))/p->DXN[IP];
                
                if(!bed)
                div += fac*(dx(IJKm1,Ip1JKm1,IJKm1,0,-1) - dx(Im1JKm1,IJKm1,Im1JKm1,-1,-1))/p->DXN[IP];
                
                if(p->j_dir==1)
                {
                div += (1.0-fac)*(dy(IJK,IJp1K,IJK,0,0) - dy(IJm1K,IJK,IJm1K,-1,0))/p->DYN[JP];
                
                if(!bed)
                div += fac*(dy(IJKm1,IJp1Km1,IJKm1,0,-1) - dy(IJm1Km1,IJKm1,IJm1Km1,-1,-1))/p->DYN[JP];
                }
                
                div += p->sigz[IJ]*(dz(IJK) - (bed ? 0.0 : dz(IJKm1)))/p->DZP[KM1];
                
                d->rhsvec.V[n] -= div/adt;
            }
        ++n;
        }
    }
    
    // continuity flux of the layers at the x and y faces near the body (after the HLL flux): the face velocity of
    // the projection on the open part, the body on the closed part
    inline void continuity_flux(lexer *p, fdm_nhf *d)
    {
        int i,j,k;
        const int marge = increment::marge;
        
        ULOOP
        if(d->FBNEAR[IJ]==1 || d->FBNEAR[Ip1J]==1)
        d->FEx[IJK] = 0.5*(d->Ds(i,j) + d->Dn(i,j))*(d->MBX[IJK]*(0.5*(d->U[IJK] + d->U[Ip1JK]) + rcx(p,d,i,j,k)) + (1.0 - d->MBX[IJK])*d->FBUX[IJK]);
        
        if(p->j_dir==1)
        VLOOP
        if(d->FBNEAR[IJ]==1 || d->FBNEAR[IJp1]==1)
        d->FEy[IJK] = 0.5*(d->De(i,j) + d->Dw(i,j))*(d->MBY[IJK]*(0.5*(d->V[IJK] + d->V[IJp1K]) + rcy(p,d,i,j,k)) + (1.0 - d->MBY[IJK])*d->FBUY[IJK]);
    }
}

#endif
