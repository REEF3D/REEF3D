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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

    /**
     * @brief Moves the particles with the flow field.
     *
     * This function is responsible for moving the particle with the flow field.
     * It does so by calculating the experienced accelerations and the resulting new velocity and position of the particle.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param pgc A reference to the ghostcell object.
     * @param PP A reference to the particles_obj object.
     */
     
void partres::advec_pic(lexer *p, fdm &a, particles_obj &PP, size_t n, sediment_fdm &s, turbulence &pturb, 
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &du, double &dv, double &dw, double alpha)
{
    double RKu,RKv,RKw;
    double uf=0.0,vf=0.0,wf=0.0;
    double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
    double ws=0;
    double topoDist=0;
       
    double dPx=0.0, dPy=0.0, dPz=0.0;
    double dTx=0.0, dTy=0.0, dTz=0.0;

    double Bx=(1.0-drho)*p->W20;
    double By=(1.0-drho)*p->W21;
    double Bz=(1.0-drho)*p->W22;
    
    double Urel=0.0,Vrel=0.0,Wrel=0.0;
    double Ts,T;
    
    double Dpx=0;
    double Dpy=0;
    double Dpz=0;
    double Uabs_rel=0.0;

    // find cell IJK
    i=p->posc_i(PX[n]);
    j=p->posc_j(PY[n]);
    k=p->posc_k(PZ[n]);

    // theta calc
    Ts = PI*pow(PP.d50,3.0)*(cellSum[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
    
        
    dTx = (stressTensor[Ip1JK] - stressTensor[Im1JK])/(p->DXP[IM1]+p->DXP[IP]);
    dTy = (stressTensor[IJp1K] - stressTensor[IJm1K])/(p->DYP[JM1]+p->DYP[JP]);
    dTz = (stressTensor[IJKp1] - stressTensor[IJKm1])/(p->DZP[KM1]+p->DZP[KP]);
    
    dPx = (a.press(i+1,j,k) - a.press(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
    dPy = (a.press(i,j+1,k) - a.press(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
    dPz = (a.press(i,j,k+1) - a.press(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
    
    //dPx = ((a.press(i+1,j,k)-a.phi(i+1,j,k)*a.ro(i+1,j,k)*fabs(p->W22)) - ((a.press(i-1,j,k)-a.phi(i-1,j,k)*a.ro(i-1,j,k)*fabs(p->W22))))/(p->DXP[IM1]+p->DXP[IP]);
    //dPy = ((a.press(i,j+1,k)-a.phi(i,j+1,k)*a.ro(i,j+1,k)*fabs(p->W22)) - ((a.press(i,j-1,k)-a.phi(i,j-1,k)*a.ro(i,j-1,k)*fabs(p->W22))))/(p->DYP[JM1]+p->DYP[JP]);
    dPz = ((a.press(i,j,k+1)-a.phi(i,j,k+1)*a.ro(i,j,k+1)*fabs(p->W22)) - ((a.press(i,j,k-1)-a.phi(i,j,k-1)*a.ro(i,j,k-1)*fabs(p->W22))))/(p->DZP[KM1]+p->DZP[KP]);

    // velocity
    velDist=0.0;
    
    uf = p->ccipol1(a.u,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
    vf = p->ccipol2(a.v,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
    wf = p->ccipol3(a.w,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);

    // relative velocity
    Urel = uf-PU[n];
    Vrel = vf-PV[n];
    Wrel = wf-PW[n];

    PP.Uf[n]=Urel;
    PP.Vf[n]=Vrel;
    PP.Wf[n]=Wrel;
    
    // drag coefficient
    //if(p->Q202==1)
    Dpx=drag_model(p,PP.d50,Urel,Ts);
    Dpy=drag_model(p,PP.d50,Vrel,Ts);
    Dpz=drag_model(p,PP.d50,Wrel,Ts);
    
    //Dp = 0.5;
        
    /*if(p->Q202==2)
    {
    Uabs_rel=sqrt(Urel*Urel + Vrel*Vrel);
    Re_p = Uabs_rel*PP.d50/(p->W2/p->W1);
    Dp = drag_coefficient(Re_p);
    }*/
    
    PP.drag[n]=Dp;

// particle force

    
    
    /*Fs = p->Q30*(p->S22-p->W1)*fabs(p->W22)/p->S22;
    

    if(Fs>du)
    du=0.0;
    
    if(Fs>dv)
    dv=0.0;
    
    cout<<"Fs: "<<Fs<<" du: "<<du<<" dv: "<<dv<<endl;*/
    
    du = Dpx*Urel + Bx - dPx/p->S22 - dTx/((Ts>1.0e10?Ts:1.0e10)*p->S22);
    dv = Dpy*Vrel + By - dPy/p->S22 - dTy/((Ts>1.0e10?Ts:1.0e10)*p->S22);
    dw = Dpz*Wrel + Bz*0.0 - dPz/p->S22*0.0 - dTz/((Ts>1.0e10?Ts:1.0e10)*p->S22);


    //cout<<"dw: "<<dw<<" Bz: "<<dTz<<" dPz: "<<-dPz/p->S22<<" dTz: "<<-dTz/((Ts>1.0e10?Ts:1.0e10)*p->S22)<<endl;

    
    // solid forcing
    double fx,fy,fz;
    if(p->S10==2)
    {
    fx = p->ccipol1c(a.fbh1,PX[n],PY[n],PZ[n])*(0.0-PU[n])/(alpha*p->dtsed); 
    fy = p->ccipol2c(a.fbh2,PX[n],PY[n],PZ[n])*(0.0-PV[n])/(alpha*p->dtsed); 
    fz = p->ccipol3c(a.fbh3,PX[n],PY[n],PZ[n])*(0.0-PW[n])/(alpha*p->dtsed); 
    
    du += fx;
    dv += fy;
    dw += fz;
    }
    
    // relax
    du *= rf(p,PX[n],PY[n]);
    dv *= rf(p,PX[n],PY[n]);
    dw *= rf(p,PX[n],PY[n]);
    

    Umax = MAX(Umax,sqrt(PU[n]*PU[n] + PV[n]*PV[n]));


    // error call
    if(PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
    cout<<"NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<Dp<<"\nTs: "<<Ts<<endl;
    cout<<"du: "<<du<<" dv: "<<dv<<" dw: "<<dw<<endl;
    cout<<"dTx: "<<dTx<<" dTy: "<<dTy<<" dTz: "<<dTz<<endl;
    exit(1);
    }
}
