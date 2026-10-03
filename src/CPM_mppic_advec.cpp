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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"part.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

// MP-PIC parcel acceleration
//   dUp/dt = Dp (Uf - Up) - grad(p)/rho_p + g - grad(Ps)/(theta rho_p)
// returns the explicit part F,G,H = Dp Uf - grad(p)/rho_p + g - grad(Ps)/(theta rho_p) (+ solid forcing)
// and the implicit drag coefficient Dpx (=Dpy=Dpz), velocity update: Up = (Up + dt F)/(1 + dt Dp)
void CPM::advec_mppic(lexer *p, fdm *a, part &P, sediment_fdm *s, turbulence *pturb,
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &F, double &G, double &H, double dt)
{
    // fluid pressure gradient
    dPx_val = p->ccipol4a(dPx,PX[n],PY[n],PZ[n]);
    dPy_val = p->j_dir==1 ? p->ccipol4a(dPy,PX[n],PY[n],PZ[n]) : 0.0;
    dPz_val = p->ccipol4a(dPz,PX[n],PY[n],PZ[n]);

    // gravity
    Bx = p->W20;
    By = p->W21;
    Bz = p->W22;

    // fluid velocity
    uf = p->ccipol1c(a->u,PX[n],PY[n],PZ[n]);
    vf = p->j_dir==1 ? p->ccipol2c(a->v,PX[n],PY[n],PZ[n]) : 0.0;
    wf = p->ccipol3c(a->w,PX[n],PY[n],PZ[n]);
    
    // S 10 1: the bed is a solid boundary for the fluid
    if(p->S10!=2)
    nearbed_velocity(p,a,PX[n],PY[n],PZ[n]);

    P.Uf[n] = uf;
    P.Vf[n] = vf;
    P.Wf[n] = wf;
    
    // relative velocity
    Urel = uf-PU[n];
    Vrel = vf-PV[n];
    Wrel = wf-PW[n];
    
    Tsval = p->ccipol4a(Ts,PX[n],PY[n],PZ[n]);
    
    // drag coefficient
    Dpx = drag_model(p,P.D[n],P.RO[n],sqrt(Urel*Urel + Vrel*Vrel + Wrel*Wrel),Tsval);
    Dpy = Dpz = Dpx;

    // explicit forces
    F = Dpx*uf - dPx_val/P.RO[n] + Bx;
    G = Dpy*vf - dPy_val/P.RO[n] + By;
    H = Dpz*wf - dPz_val/P.RO[n] + Bz;
    
    // inter-particle normal stress
    if(p->Q12>0)
    {
    dTx_val = p->ccipol4a(dTx,PX[n],PY[n],PZ[n]);
    dTy_val = p->j_dir==1 ? p->ccipol4a(dTy,PX[n],PY[n],PZ[n]) : 0.0;
    dTz_val = p->ccipol4a(dTz,PX[n],PY[n],PZ[n]);
    
    double Tdiv = P.RO[n]*MAX(Tsval, 0.5*theta_bed);
    
    F -= dTx_val/Tdiv;
    G -= dTy_val/Tdiv;
    H -= dTz_val/Tdiv;
    }
    
    if(p->j_dir==0)
    G=0.0;
    
    P.Test[n] = Tsval;

    // solid forcing
    double fx,fy,fz;
    fx = p->ccipol1c(a->fbh1,PX[n],PY[n],PZ[n])*(0.0-PU[n])/dt;
    fy = p->ccipol2c(a->fbh2,PX[n],PY[n],PZ[n])*(0.0-PV[n])/dt;
    fz = p->ccipol3c(a->fbh3,PX[n],PY[n],PZ[n])*(0.0-PW[n])/dt;

    F += fx;
    G += fy;
    H += fz;

    // relax
    if(p->Q73>0)
    {
    double r = rf(p,PX[n],PY[n]);
    F *= r;
    G *= r;
    H *= r;
    }

    // error call
    if(F!=F || G!=G || H!=H || PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
        cout<<"CPM NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<Dpx<<"\nTs: "<<Tsval<<endl;
        cout<<"F: "<<F<<" G: "<<G<<" H: "<<H<<endl;
        cout<<"dTx: "<<dTx_val<<" dTy: "<<dTy_val<<" dTz: "<<dTz_val<<endl;
        exit(1);
    }
}

/*--------------------------------------------------------------------
fluid velocity at parcels inside the bed for S 10 1, where the bed (topo<0) is a
solid boundary for the fluid and the bed surface lies on the grid scale:

  - surface layer, depth below the bed surface delta < 0.5 h: exposed grains,
    fluid velocity sampled 0.5 h above the bed surface along the bed normal and
    scaled to the grain with the logarithmic law of the wall,
        u_g = u_ref ln(30 d50/ks)/ln(30 (0.5h)/ks),   ks = S 21 d50
  - deeper: the pore water moves with the grains (local solid velocity)
  - linear transition between 0.5 h and h
--------------------------------------------------------------------*/
void CPM::nearbed_velocity(lexer *p, fdm *a, double xp, double yp, double zp)
{
    double topo = p->ccipol4_b(a->topo,xp,yp,zp);
    
    if(topo>=0.0)
    return;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
    
    double h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
    double delta = -topo;
    
    // pore water with the grains
    double ug = p->ccipol4a(Us,xp,yp,zp);
    double vg = p->j_dir==1 ? p->ccipol4a(Vs,xp,yp,zp) : 0.0;
    double wg = p->ccipol4a(Ws,xp,yp,zp);
    
    double w = 0.0;
    
    if(delta<0.5*h)
    w = 1.0;
    
    else if(delta<h)
    w = (h-delta)/(0.5*h);
    
    if(w>0.0)
    {
        // bed normal from the topo level set
        double nx = (p->ccipol4_b(a->topo,xp+0.5*h,yp,zp) - p->ccipol4_b(a->topo,xp-0.5*h,yp,zp));
        double ny = p->j_dir==1 ? (p->ccipol4_b(a->topo,xp,yp+0.5*h,zp) - p->ccipol4_b(a->topo,xp,yp-0.5*h,zp)) : 0.0;
        double nz = (p->ccipol4_b(a->topo,xp,yp,zp+0.5*h) - p->ccipol4_b(a->topo,xp,yp,zp-0.5*h));
        double nm = sqrt(nx*nx + ny*ny + nz*nz);
        
        if(nm>1.0e-10)
        {
            nx/=nm;
            ny/=nm;
            nz/=nm;
        }
        else
        {
            nx=ny=0.0;
            nz=1.0;
        }
        
        double zref = 0.5*h;
        double xr = xp + nx*(delta+zref);
        double yr = yp + ny*(delta+zref);
        double zr = zp + nz*(delta+zref);
        
        double ur = p->ccipol1c(a->u,xr,yr,zr);
        double vr = p->j_dir==1 ? p->ccipol2c(a->v,xr,yr,zr) : 0.0;
        double wr = p->ccipol3c(a->w,xr,yr,zr);
        
        double ks = MAX(p->S21*p->S20, 1.0e-6);
        double fac = log(MAX(30.0*p->S20/ks,1.0+1.0e-6))/log(MAX(30.0*zref/ks,1.0+1.0e-6));
        fac = MAX(0.0,MIN(fac,1.0));
        
        uf = w*fac*ur + (1.0-w)*ug;
        vf = w*fac*vr + (1.0-w)*vg;
        wf = w*fac*wr + (1.0-w)*wg;
    }
    else
    {
        uf = ug;
        vf = vg;
        wf = wg;
    }
}
