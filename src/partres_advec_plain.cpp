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
#include"part.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

void partres::advec_plain(lexer *p, fdm *a, part &P, sediment_fdm *s, turbulence *pturb, 
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &F, double &G, double &H, double alpha)
{
    // find cell IJK
    i=p->posc_i(PX[n]);
    j=p->posc_j(PY[n]);
    k=p->posc_k(PZ[n]);
   
    // velocity
    uf = p->ccipol1(a->u,PX[n],PY[n],PZ[n]);
    vf = p->ccipol2(a->v,PX[n],PY[n],PZ[n]);
    wf = p->ccipol3(a->w,PX[n],PY[n],PZ[n]);
    
    // relative velocity
    Urel = uf-PU[n];
    Vrel = vf-PV[n];
    Wrel = wf-PW[n];
    
    Uabs_rel=sqrt(Urel*Urel + Vrel*Vrel);

    P.Uf[n]=uf;
    P.Vf[n]=vf;
    P.Wf[n]=wf;
    
    DragCoeff = 0.5;
    
    // particle force
    // acceleration
    Fd = p->W1 * DragCoeff * PI/8.0 * pow(P.d50,2)  * pow(Uabs_rel,2.0);
    
    // relax
    Fd *= rf(p,PX[n],PY[n]);
    
    // resistance force
    Fs = p->Q30*(p->S22-p->W1)*fabs(p->W22)*PI*pow(P.d50, 3.0)/6.0;
     
    F_tot = Fd-Fs;//*s.reduce(i,j);
    
    F_tot = MAX(F_tot,0.0);
    
    //cout<<"Fd: "<<Fd<<" Fs: "<<Fs<<" F_tot: "<<F_tot<<" "<<P.d50<<" "<<DragCoeff<<endl;


// particle force
    F = F_tot /(p->S22*PI*pow(P.d50,3.0)/6.0) * Urel/(fabs(Uabs_rel)>1.0e-10?fabs(Uabs_rel):1.0e20);
    G = F_tot /(p->S22*PI*pow(P.d50,3.0)/12.0) * Vrel/(fabs(Uabs_rel)>1.0e-10?fabs(Uabs_rel):1.0e20);
    H = 0.0;
    
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
    
    if(PX[n]<1.9)
    F=G=H=0.0;
    
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
