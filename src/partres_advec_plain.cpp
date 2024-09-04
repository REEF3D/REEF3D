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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "partres.h"
#include "particles_obj.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

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

    
    double RKtimeStep = 0.5*p->dtsed; // To be modified with alpha?


    double Du=0,Dv=0,Dw=0;
    double thetas=0;
    double DragCoeff=0;
    double dU=0;

    i=p->posc_i(PX[n]);
    j=p->posc_j(PY[n]);
    k=p->posc_k(PZ[n]);

    thetas=theta_s(p,a,PP,i,j,k);
        

    topoDist=p->ccipol4(a.topo,PX[n],PY[n],PZ[n]);

    if(topoDist<velDist*p->DZP[KP])
    {
        u=p->ccipol1c(a.u,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]-topoDist);
        v=p->ccipol2c(a.v,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]-topoDist);
        // w=p->ccipol3c(a.w,PX[n],PY[n],PZ[n]+velDist*p->DZP[KP]-topoDist);
        if(debugPrint)
        {
                cout<<PZ[n]+velDist*p->DZP[KP]-topoDist<<endl;
                debugPrint=false;
        }
    }
    
    else
    {
        u=p->ccipol1c(a.u,PX[n],PY[n],PZ[n]);
        v=p->ccipol2c(a.v,PX[n],PY[n],PZ[n]);
        // w=p->ccipol3c(a.w,PX[n],PY[n],PZ[n]);
    }

    // PP.Uf[n]=u;
    // PP.Vf[n]=v;
    // PP.Wf[n]=w;
        
    Du=u-PU[n];
    Dv=v-PV[n];
    // Dw=w-PW[n];

    switch (p->Q202)
    {
        case 0:
        {
            DragCoeff=drag_model(p,PP.d50,Du,Dv,Dw,thetas);
            break;
        }
        
        case 1:
        {
            relative_velocity(p,a,PP,n,Du,Dv,Dw);
            dU=sqrt(Du*Du+Dv*Dv+Dw*Dw);
            const double Re_p = dU*PP.d50/(p->W2/p->W1);
            DragCoeff=drag_coefficient(Re_p);
            break;
        }
    }
    PP.drag[n]=DragCoeff;

    // Acceleration
    switch (p->Q202)
    {
        case 0:
        {
            du=DragCoeff*Du;
            dv=DragCoeff*Dv;
            // dw=DragCoeff*Dw;
            du+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22);
            dv+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22);
            // dw+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22);
            break;
        }
        
        case 1:
        {
            const double Fd = DragCoeff * PI/8.0 * pow(PP.d50,2) * p->W1 * pow(dU,2);
            DragCoeff = Fd /(p->S22*PI*pow(PP.d50,3.0)/6.0);
            du=DragCoeff;
            dv=DragCoeff;
            // dw=DragCoeff;
                break;
        }
    }


    // solid forcing
    double fx,fy,fz;
    
    fx = p->ccipol1c(a.fbh1,PX[n],PY[n],PZ[n])*(0.0-PU[n])/(alpha*p->dtsed); 
    fy = p->ccipol2c(a.fbh2,PX[n],PY[n],PZ[n])*(0.0-PV[n])/(alpha*p->dtsed); 
    fz = p->ccipol3c(a.fbh3,PX[n],PY[n],PZ[n])*(0.0-PW[n])/(alpha*p->dtsed); 
    
    du += fx;
    dv += fx;
    dw += fx;
    
    if(PU[n]!=PU[n] || PV[n]!=PV[n] || PW[n]!=PW[n])
    {
    cout<<"NaN detected.\nDu: "<<Du<<" Dv: "<<Dv<<" Dw: "<<Dw<<"\nDrag: "<<DragCoeff<<endl;
    exit(1);
    }
}
