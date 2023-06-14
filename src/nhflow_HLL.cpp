/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_HLL.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"
#include"nhflow_flux_face_cds2.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_reconstruct_WENO.h"

nhflow_HLL::nhflow_HLL (lexer *p, patchBC_interface *ppBC)
{
    pBC = ppBC;
    
    pflux = new nhflow_flux_face_cds2(p);
    
    if(p->A543==2)
    precon = new nhflow_reconstruct_hires(p,ppBC);
    
    if(p->A543==4)
    precon = new nhflow_reconstruct_weno(p,ppBC);
}

nhflow_HLL::~nhflow_HLL()
{
}

void nhflow_HLL::start(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W)
{
    // reconstruct eta
    precon->reconstruct_2D(p, pgc, d, eta, ETAs, ETAn, ETAe, ETAw);
    
   
    
        if(ipol==1)
        LOOP
        d->F[IJK]+=aij(p,d,F,1,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==2)
        LOOP
        d->G[IJK]+=aij(p,d,F,2,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==3)
        LOOP
        d->H[IJK]+=aij(p,d,F,3,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==4)
        LOOP
        d->L[IJK]+=aij(p,d,F,4,U,V,W,p->DXN,p->DYN,p->DZN);
}

double nhflow_HLL::aij_U(lexer* p,fdm_nhf* d, double *F, int ipol, double *UVEL, double *VVEL, double *WVEL, double *DX,double *DY, double *DZ)
{
    // reconstruct
    /*
    EtaxL
    EtaxR
    
    UxL
    UxR
    DUxL
    DUxR
    
    UyL
    UyR
    DUyL
    DUyR
    
    UzL
    UzR
    DUxL
    DUxR
    
    ----
    
    FxL
    FxR
    FyL
    FyR
    Fz  
    
    
    FxL(i,j,k) = DUxL(i,j,k)*UxL(i,j,k)+0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx0(i,j))
    FxR(i,j,k) = DUxR(i,j,k)*UxR(i,j,k)+0.5*Grav*(EtaxR(i,j)*EtaxR(i,j)+2.0*EtaxR(i,j)*hfx0(i,j))
    
    FyL(i,j,k) = DyL(i,j)*UyL(i,j,k)*VyL(i,j,k)
    FyR(i,j,k) = DyR(i,j)*UyR(i,j,k)*VyR(i,j,k)
    
    Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
    */
    
    
    double Ss,Sn,Se,Sw;
    double USx,USy;
    double Ds,Dn,De,Dw;
    double DSx,DSy;
    double denom;
    
    LOOP
    DV[IJK] = UVEL[IJK]*d->WL[IJ];
    
     // reconstruct U and V
    precon->reconstruct_3D(p, pgc, d, DV, DV, Fs, Fn, Fe, Fw);
    precon->reconstruct_3D(p, pgc, d, UVEL, UVEL, Us, Un, Ve, Vw);
    
    // velocity depth
    LOOP
    {
    DU[IJK] = U[IJK]*d->WL[IJ];
    DV[IJK] = V[IJK]*d->WL[IJ];
    }
    
    pgc->start1V(p,DU,10);
    pgc->start2V(p,DV,11);
         
    
    
    // HLL flux
    ULOOP
    {
    // water level       
    Ds = ETAs(i,j) + p->wd - d->bed(i,j);
    Dn = ETAn(i,j) + p->wd - d->bed(i,j);
    
    Ds = MAX(0.00005, Ds);
    Dn = MAX(0.00005, Dn);
    
    // Us
    USx = 0.5*(Us[IJK]+Un[IJK]) + sqrt(9.81*Ds) - sqrt(9.81*Dn);
    DSx = 0.5*(sqrt(9.81*Ds) + sqrt(9.81*Dn)) + 0.25*(Us[IJK] - Un[IJK]);
    
    // wave speed
    Ss = MIN(Us[IJK] - sqrt(9.81*Ds), USx - DSx);
    Sn = MAX(Un[IJK] + sqrt(9.81*Dn), USx + DSx);
    
    if(p->wet[Ip1J]==0)
    {
    Ss = Us[IJK] - sqrt(9.81*Ds);
    Sn = Us[IJK] + 2.0*sqrt(9.81*Ds);
    }
    
    if(p->wet[Im1J]==0)
    {
    Ss = Un[IJK] - 2.0*sqrt(9.81*Dn);
    Sn = Un[IJK] + sqrt(9.81*Dn);
    }
    
    if(p->wet[IJ]==0)
    {
    Ss=Sn=0.0;
    USx=0.0;
    }
    
        // final flux x-dir
        if(Ss>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(Sn<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        {
        denom = Sn-Ss;
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fx[IJK] = (Sn*Fs[IJK] - Ss*Fn[IJK] + Sn*Ss*(Dn - Ds))/denom;
        }
    }
    
    
    VLOOP
    {
    // water level       
    De = ETAe(i,j)  + p->wd - d->bed(i,j);
    Dw = ETAw(i,j)  + p->wd - d->bed(i,j);
    
    De = MAX(0.00005, De);
    Dw = MAX(0.00005, Dw);
    
    // Us
    USy = 0.5*(Ve[IJK]+Vw[IJK]) + sqrt(9.81*De) - sqrt(9.81*Dw);
    DSy = 0.5*(sqrt(9.81*De) + sqrt(9.81*Dw)) + 0.25*(Ve[IJK] - Vw[IJK]);
    
    // wave speed
    Se = MIN(Ve[IJK] - sqrt(9.81*De), USy - DSy);
    Sw = MAX(Vw[IJK] + sqrt(9.81*Dw), USy + DSy);
    
    if(p->wet[IJp1]==0)
    {
    Se = Ve[IJK] - sqrt(9.81*De);
    Sw = Ve[IJK] + 2.0*sqrt(9.81*De);
    }
    
    if(p->wet[IJm1]==0)
    {
    Se = Vw[IJK] - 2.0*sqrt(9.81*Dw);
    Sw = Vw[IJK] + sqrt(9.81*Dw);
    }
    
    if(p->wet[IJ]==0)
    {
    Se=Sw=0.0;
    USy=0.0;
    }

        // final flux y-dir
        if(Se>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(Sw<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        {
        denom = Sw-Se;
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fy[IJK] = (Sw*Fe[IJK] - Se*Fw[IJK] + Sw*Se*(Dw - De))/denom;
        }
    }

    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,10);
	
    
}



double nhflow_HLL::aij(lexer* p,fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, double *DX,double *DY, double *DZ)
{
    udir=vdir=wdir=0.0;
    
        // convective flux
        pflux->u_flux(d,ipol,U,ivel1,ivel2);
        pflux->v_flux(d,ipol,V,jvel1,jvel2);
        pflux->w_flux(d,ipol,d->omegaF,kvel1,kvel2);

    
    // x-dir
    if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;
    
    dx =     udir*(ivel2*F[IJK] - ivel1*F[Im1JK])/DX[IM1] 
    
    +   (1.0-udir)*(ivel2*F[Ip1JK] - ivel1*F[IJK])/DX[IP]; 
    
    
    // y-dir
    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;
    
    dy =     vdir*(jvel2*F[IJK] - jvel1*F[IJm1K])/DY[JM1] 
    
    +   (1.0-vdir)*(jvel2*F[IJp1K] - jvel1*F[IJK])/DY[JP]; 
    
    
    // z-dir
    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;
    
    dz =     wdir*(kvel2*F[IJK] - kvel1*F[IJKm1])/DZ[KM1] 
    
    +   (1.0-wdir)*(kvel2*F[IJKp1] - kvel1*F[IJK])/DZ[KP]; 
    
    
    L = -dx-dy-dz;
    
    return L;
}

