/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"part.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

void partres::advec_mppic(lexer *p, fdm *a, part &P, sediment_fdm *s, turbulence *pturb, 
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &F, double &G, double &H, double alpha)
{
    // find cell IJK
    i=p->posc_i(PX[n]);
    j=p->posc_j(PY[n]);
    k=p->posc_k(PZ[n]);
    
    // pressure gradient
    dPx_val = p->ccipol4a(dPx,PX[n],PY[n],PZ[n]);
    dPy_val = p->ccipol4a(dPy,PX[n],PY[n],PZ[n]);
    dPz_val = p->ccipol4a(dPz,PX[n],PY[n],PZ[n]);
   
    // buouancy
    Bx = p->W20;
    By = p->W21;
    Bz = p->W22;
    
    // inter-particle stress
    Tsval = p->ccipol4a(Ts,PX[n],PY[n],PZ[n]);
    
    dTx_val = p->ccipol4a(dTx,PX[n],PY[n],PZ[n]);
    dTy_val = p->ccipol4a(dTy,PX[n],PY[n],PZ[n]);
    dTz_val = p->ccipol4a(dTz,PX[n],PY[n],PZ[n]);
    
    // velocity
    uf = p->ccipol1(a->u,PX[n],PY[n],PZ[n]);
    vf = p->ccipol2(a->v,PX[n],PY[n],PZ[n]);
    wf = p->ccipol3(a->w,PX[n],PY[n],PZ[n]);
    
    // relative velocity
    Urel = uf-PU[n];
    Vrel = vf-PV[n];
    Wrel = wf-PW[n];

    P.Uf[n]=uf;
    P.Vf[n]=vf;
    P.Wf[n]=wf;
    
    // drag coefficient
    Dpx=drag_model(p,P.D[n],P.RO[n],Urel,Tsval);
    Dpy=drag_model(p,P.D[n],P.RO[n],Vrel,Tsval);
    Dpz=drag_model(p,P.D[n],P.RO[n],Wrel,Tsval);
    
// particle force
    F = 0.0*Dpx*Urel - 0.0*dPx_val/P.RO[n] + Bx - 0.0*dTx_val/P.RO[n];
    G = 0.0*Dpy*Vrel - 0.0*dPy_val/P.RO[n] + By - 0.0*dTy_val/P.RO[n];
    H = 0.0*Dpz*Wrel - dPz_val/P.RO[n] + Bz - dTz_val/P.RO[n];
    
    // solid forcing
    double fx,fy,fz;
    if(p->S10==2)
    {
    fx = p->ccipol1c(a->fbh1,PX[n],PY[n],PZ[n])*(0.0-PU[n])/(alpha*p->dtsed); 
    fy = p->ccipol2c(a->fbh2,PX[n],PY[n],PZ[n])*(0.0-PV[n])/(alpha*p->dtsed); 
    fz = p->ccipol3c(a->fbh3,PX[n],PY[n],PZ[n])*(0.0-PW[n])/(alpha*p->dtsed); 
    
    F += fx;
    G += fy;
    H += fz;
    }
    
    // relax
    F *= rf(p,PX[n],PY[n]);
    G *= rf(p,PX[n],PY[n]);
    H *= rf(p,PX[n],PY[n]);
    
    //if(PX[n]<1.9)
    //F=G=H=0.0;
    
    /*F=0.5;
    G=0.0;
    H=0.0;*/
    
    P.Test[n] = dTz_val/P.RO[n];
    
    Umax = MAX(Umax,sqrt(PU[n]*PU[n] + PV[n]*PV[n]));
    
    // error call
    if(PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
    cout<<"NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<Dp<<"\nTs: "<<Tsval<<endl;
    cout<<"F: "<<F<<" G: "<<G<<" H: "<<H<<endl;
    cout<<"dTx: "<<dTx_val<<" dTy: "<<dTy_val<<" dTz: "<<dTz_val<<endl;
    exit(1);
    }    
}


