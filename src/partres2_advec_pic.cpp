/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"partres2.h"
#include"part.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

void partres2::advec_pic(lexer *p, fdm *a, part &P, sediment_fdm *s, turbulence *pturb, 
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &F, double &G, double &H, double alpha)
{
    
    Bx=(1.0-p->W1/p->S22)*p->W20;
    By=(1.0-p->W1/p->S22)*p->W21;
    Bz=(1.0-p->W1/p->S22)*p->W22;
    

    // find cell IJK
    i=p->posc_i(PX[n]);
    j=p->posc_j(PY[n]);
    k=p->posc_k(PZ[n]);
    
    //cout<<n<<" PX[n]: "<<PX[n]<<" PY[n]: "<<PY[n]<<" PZ[n]: "<<PZ[n]<<endl;
    
    //cout<<n<<" P.index: "<<P.index<<" i: "<<i<<" j: "<<j<<" k: "<<k<<endl;

    // theta calc
    dTx = (Tau(i+1,j,k) - Tau(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
    dTy = (Tau(i,j+1,k) - Tau(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
    dTz = (Tau(i,j,k+1) - Tau(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
    
    //cout<<"dTx: "<<dTx<<" dTy: "<<dTy<<" dTz: "<<dTz<<endl;
    
    dPx = (a->press(i+1,j,k) - a->press(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
    dPy = (a->press(i,j+1,k) - a->press(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
    dPz = (a->press(i,j,k+1) - a->press(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
    
    //dPx = ((a->press(i+1,j,k)-a->phi(i+1,j,k)*a->ro(i+1,j,k)*fabs(p->W22)) - ((a->press(i-1,j,k)-a->phi(i-1,j,k)*a->ro(i-1,j,k)*fabs(p->W22))))/(p->DXP[IM1]+p->DXP[IP]);
    //dPy = ((a->press(i,j+1,k)-a->phi(i,j+1,k)*a->ro(i,j+1,k)*fabs(p->W22)) - ((a->press(i,j-1,k)-a->phi(i,j-1,k)*a->ro(i,j-1,k)*fabs(p->W22))))/(p->DYP[JM1]+p->DYP[JP]);
    dPz = ((a->press(i,j,k+1)-a->phi(i,j,k+1)*a->ro(i,j,k+1)*fabs(p->W22)) - ((a->press(i,j,k-1)-a->phi(i,j,k-1)*a->ro(i,j,k-1)*fabs(p->W22))))/(p->DZP[KM1]+p->DZP[KP]);

    // velocity
    velDist=0.0;
    
    uf = p->ccipol1(a->u,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
    vf = p->ccipol2(a->v,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
    wf = p->ccipol3(a->w,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
    
    //cout<<"uf: "<<uf<<" vf: "<<vf<<" wf: "<<wf<<endl;

    // relative velocity
    Urel = uf-PU[n];
    Vrel = vf-PV[n];
    Wrel = wf-PW[n];

    P.Uf[n]=uf;
    P.Vf[n]=vf;
    P.Wf[n]=wf;
    
    // drag coefficient
    Dpx=drag_model(p,P.D[n],P.RO[n],Urel,Ts);
    Dpy=drag_model(p,P.D[n],P.RO[n],Vrel,Ts);
    Dpz=drag_model(p,P.D[n],P.RO[n],Wrel,Ts);
    
    //cout<<"Dpx: "<<Dpx<<" Dpy: "<<Dpy<<" Dpz: "<<Dpz<<endl;
    

// particle force
    
    F = Dpx*Urel + Bx - dPx/p->S22 - dTx/((Ts>1.0e10?Ts:1.0e10)*p->S22);
    G = Dpy*Vrel + By - dPy/p->S22 - dTy/((Ts>1.0e10?Ts:1.0e10)*p->S22);
    H = Dpz*Wrel + Bz*0.0 - dPz/p->S22*0.0 - dTz/((Ts>1.0e10?Ts:1.0e10)*p->S22);

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
    

    Umax = MAX(Umax,sqrt(PU[n]*PU[n] + PV[n]*PV[n]));
    
    
    //cout<<"F: "<<F<<" G: "<<G<<" H: "<<H<<endl;
    
    /*F=0.5;
    G=0.0;
    H=0.0;*/

    // error call
    if(PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
    cout<<"NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<Dp<<"\nTs: "<<Ts<<endl;
    cout<<"F: "<<F<<" G: "<<G<<" H: "<<H<<endl;
    cout<<"dTx: "<<dTx<<" dTy: "<<dTy<<" dTz: "<<dTz<<endl;
    exit(1);
    }
    
    
    
    
}