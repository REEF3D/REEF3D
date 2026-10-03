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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans_definitions.h"

// Porous floating body (X 16 n d50 alpha beta)
//
// Instead of forcing the fluid inside the body to the rigid-body velocity, the body acts
// as a moving porous medium: a Darcy-Forchheimer resistance on the relative (Darcy) velocity
//
//      f = -H * K * (u - u_fb),     K = Apor*visc + Bpor*|u - u_fb|
//
// is applied, integrated implicitly over the substep so that large K is stable:
//
//      u^{n+1} - u_fb = (u* - u_fb) / (1 + a dt CPOR K)
//      FX = -H K/(1 + a dt CPOR K) (u* - u_fb)
//
// K -> inf recovers the rigid direct forcing, K -> 0 a transparent body.
// The reaction -rho*FX*dV is accumulated as the drag on the skeleton (Xd..Nd) and added to
// the pressure load in force_calc_stl(), where the envelope pressure is scaled by (1-n).
//
// The drag couples fluid and body stiffly (rate ~ rho*K/rho_bulk). The fluid side is implicit
// through Kimp; the body side is made linearly implicit in porous_damping_nhflow() with the
// linearised drag coefficients Dpor_t, Dpor_r collected here.

void sixdof_obj_nhflow::update_forcing_nhflow_porous(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, int iter)
{
    double H,uf,vf,wf;
    double du,dv,dw,urel;
    double Kval,Kimp,fx,fy,fz;
    double dV,rx,ry,rz;
    
    const double dtsub = alpha[iter]*p->dt;
    
    // porosity field follows the body (also stamps d->FHB)
    porosity_nhflow(p,d,pgc);
    
    Xd=Yd=Zd=Kd=Md=Nd=0.0;
    Dpor_t=Dpor_r[0]=Dpor_r[1]=Dpor_r[2]=0.0;
    
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        if(H<1.0e-12)
        continue;
        
        uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
        vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
        wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
        
        // relative Darcy velocity
        du = U[IJK] - uf;
        dv = V[IJK] - vf;
        dw = W[IJK] - wf;
        
        if(p->j_dir==0)
        dv = 0.0;
        
        urel = sqrt(du*du + dv*dv + dw*dw);
        
        Kval = Apor_fb*d->VISC[IJK] + Bpor_fb*urel;
        
        Kimp = Kval/(1.0 + dtsub*CPORNH*Kval);
        
        fx = -H*Kimp*du;
        fy = -H*Kimp*dv;
        fz = -H*Kimp*dw;
        
        FX[IJK] += fx;
        FY[IJK] += fy;
        FZ[IJK] += fz;
        
        // reaction on the body
        dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
        
        rx = p->pos_x() - c_(0);
        ry = p->pos_y() - c_(1);
        rz = p->pos_z() - c_(2);
        
        Xd -= p->W1*fx*dV;
        Yd -= p->W1*fy*dV;
        Zd -= p->W1*fz*dV;
        
        Kd -= p->W1*(ry*fz - rz*fy)*dV;
        Md -= p->W1*(rz*fx - rx*fz)*dV;
        Nd -= p->W1*(rx*fy - ry*fx)*dV;
        
        // linearised d(drag)/d(u_fb), translation and rotation about the CoG
        Dpor_t    += p->W1*H*Kimp*dV;
        Dpor_r[0] += p->W1*H*Kimp*(ry*ry + rz*rz)*dV;
        Dpor_r[1] += p->W1*H*Kimp*(rx*rx + rz*rz)*dV;
        Dpor_r[2] += p->W1*H*Kimp*(rx*rx + ry*ry)*dV;
    }
    
    Xd = pgc->globalsum(Xd);
    Yd = pgc->globalsum(Yd);
    Zd = pgc->globalsum(Zd);
    Kd = pgc->globalsum(Kd);
    Md = pgc->globalsum(Md);
    Nd = pgc->globalsum(Nd);
    
    Dpor_t    = pgc->globalsum(Dpor_t);
    Dpor_r[0] = pgc->globalsum(Dpor_r[0]);
    Dpor_r[1] = pgc->globalsum(Dpor_r[1]);
    Dpor_r[2] = pgc->globalsum(Dpor_r[2]);
    
    if(p->j_dir==0)
    Yd = Kd = Nd = 0.0;
}

void sixdof_obj_nhflow::porosity_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // n(x,t) = 1 - H_fb (1 - n_fb)
    // without static VRANS structures (B 200 0) nhflow_forcing::reset() clears POR to 1 before
    // the bodies stamp it; with static structures (B 200 1) vrans_nhflow_f::update() resets POR at the
    // start of every stage and applies the same stamp, here we only take the minimum.
    double H;
    
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        // POR is reset to 1 in nhflow_forcing::reset() (B 200 0) or vrans_nhflow_f::update()
        // (B 200 1) -> MIN lets several porous bodies coexist
        d->POR[IJK] = MIN(d->POR[IJK], 1.0 - H*(1.0 - p->X16_n));
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0);
    }
    
    pgc->start5Vfull(p,d->POR,1);
    pgc->start5V(p,d->FHB,50);
}

void sixdof_obj_nhflow::porous_damping_nhflow(lexer *p, int iter)
{
    // Linearly implicit treatment of the drag in the body equations:
    //   (m + a dt D) dp = a dt F   ->   F_eff = F * m/(m + a dt D)
    // Steady states are unchanged, only the stiff drag relaxation is damped.
    // Rotation uses the diagonal of the inertia tensor (exact for 2D pitch, approximate for
    // large 3D rotations).
    const double dtsub = alpha[iter]*p->dt;
    
    Ffb_ *= Mass_fb/(Mass_fb + dtsub*Dpor_t);
    
    for(int q=0; q<3; ++q)
    if(I_(q,q)>1.0e-20)
    Mfb_(q) *= I_(q,q)/(I_(q,q) + dtsub*Dpor_r[q]);
}
