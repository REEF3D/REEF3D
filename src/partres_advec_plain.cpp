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
#include"sediment_fdm.h"
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
     
void partres::advec_plain(lexer *p, fdm &a, particles_obj &PP, size_t n, sediment_fdm &s, turbulence &pturb, 
                        double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW,
                        double &du, double &dv, double &dw, double alpha)
{
    double RKu,RKv,RKw;
    double u=0,v=0,w=0;
    double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
    double ws=0;
    double topoDist=0;
       
    double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
    double stressDivX=0, stressDivY=0, stressDivZ=0;
    double netBuoyX=0, netBuoyY=0, netBuoyZ=0;
    // double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;

    bool limited = true;
    bool debugPrint = false;
    bool bedLoad = false;
    bool shearVel = true;

    double Urel=0.0,Vrel=0.0,Wrel=0.0;
    double thetas=0;
    double DragCoeff=0;
    double Uabs_rel=0.0;
    double Uabs=0.0;

    i=p->posc_i(PX[n]);
    j=p->posc_j(PY[n]);
    k=p->posc_k(PZ[n]);
        
    topoDist=p->ccipol4_a(a.topo,PX[n],PY[n],PZ[n]);
    
    velDist=0.6;

    
    if(topoDist<velDist*p->DZP[KP])
    {
        u=p->ccipol1c(a.u,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]-topoDist);
        v=p->ccipol2c(a.v,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]-topoDist);
    }
    
    
    else
    {
        u=p->ccipol1c(a.u,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
        v=p->ccipol2c(a.v,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]);
        // w=p->ccipol3c(a.w,PX[n],PY[n],PZ[n]);
    }

    // relative velocity
    Urel=u-PU[n];
    Vrel=v-PV[n];


    // particle force
    if(p->Q202==1)
    {
    DragCoeff=drag_model(p,PP.d50,Urel,Vrel,Wrel,thetas);
    
    du=DragCoeff*Urel;
    dv=DragCoeff*Vrel;
    }
        
    if(p->Q202==2)
    {
    relative_velocity(p,a,PP,n,Urel,Vrel,Wrel);
    Uabs_rel=sqrt(Urel*Urel+Vrel*Vrel+Wrel*Wrel);
    Re_p = Uabs_rel*PP.d50/(p->W2/p->W1);
    DragCoeff=drag_coefficient(Re_p);
         
    // acceleration
    Fd = p->W1 * DragCoeff * PI/8.0 * pow(PP.d50,2)  * pow(Uabs_rel,2.0);
    
    Fd *= rf(p,PX[n],PY[n]);
    
    Fs = (p->S22-p->W1)*fabs(p->W22)*PI*pow(PP.d50, 3.0)*0.88/6.0;
    
    Ft = sin(s.teta(i,j))*(p->S22-p->W1)*fabs(p->W22)*PI*pow(PP.d50, 3.0);
    
    F_tot = Fd-Fs+Ft;
    
    F_tot = MAX(F_tot,0.0);
    
    //cout<<"Fd: "<<Fd<<" Fs: "<<Fs<<" F_tot: "<<F_tot<<" "<<PP.d50<<" "<<DragCoeff<<endl;

    du = F_tot /(p->S22*PI*pow(PP.d50,3.0)/6.0) * Urel/(fabs(Uabs_rel)>1.0e-10?fabs(Uabs_rel):1.0e20);
    dv = F_tot /(p->S22*PI*pow(PP.d50,3.0)/6.0) * Vrel/(fabs(Uabs_rel)>1.0e-10?fabs(Uabs_rel):1.0e20);
    
    PP.shear_eff[n]=Fd;
    PP.shear_crit[n]=Fs;
    }

    PP.drag[n]=DragCoeff;
    
    Uabs=sqrt(PU[n]*PU[n]+PV[n]*PV[n]);
    
    if(Uabs<0.1)
    PP.Flag[n]=0;
    
    if(Uabs>=0.1)
    PP.Flag[n]=1;

    // solid forcing
    double fx,fy,fz;
    /*
    fx = p->ccipol1c(a.fbh1,PX[n],PY[n],PZ[n])*(0.0-PU[n])/(alpha*p->dtsed); 
    fy = p->ccipol2c(a.fbh2,PX[n],PY[n],PZ[n])*(0.0-PV[n])/(alpha*p->dtsed); 
    fz = p->ccipol3c(a.fbh3,PX[n],PY[n],PZ[n])*(0.0-PW[n])/(alpha*p->dtsed); 
    
    du += fx;
    dv += fy;
    dw += fz;*/
    
    // relax
    du *= rf(p,PX[n],PY[n]);
    dv *= rf(p,PX[n],PY[n]);
    dw *= rf(p,PX[n],PY[n]);
    
    dw=0.0;
    
    
    if(PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
    cout<<"NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<DragCoeff<<endl;
    exit(1);
    }
}
