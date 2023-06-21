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
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_flux_face_cds2.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"

#define WLVL (fabs(d->WL_n1(i,j))>1.0e-20?d->WL_n1(i,j):1.0e20)

nhflow_HLL::nhflow_HLL (lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) : ETAs(p),ETAn(p),ETAe(p),ETAw(p),
                                                                              Ds(p),Dn(p),De(p),Dw(p)
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_face_cds2(p);
 
    precon = new nhflow_reconstruct_hires(p,ppBC);
    
    pss = new nhflow_signal_speed(p);
    
    p->Darray(Fx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fy,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fz,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Fs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fw,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Us,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Un,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ue,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Uw,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ub,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ut,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Vs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ve,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vw,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vt,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Ws,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Wn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(We,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ww,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Wb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Wt,p->imax*p->jmax*(p->kmax+2));
}

nhflow_HLL::~nhflow_HLL()
{
}

void nhflow_HLL::precalc(lexer* p, fdm_nhf* d, double *F, int ipol, double *UVEL, double *VVEL, double *WVEL, slice &eta)
{
    // reconstruct eta
    precon->reconstruct_2D(p, pgc, d, eta, ETAs, ETAn, ETAe, ETAw);
    
    SLICELOOP1
    {
    // water level       
    Ds(i,j) = ETAs(i,j) + 0.5*(d->depth(i,j) + d->depth(i-1,j));
    Dn(i,j) = ETAn(i,j) + 0.5*(d->depth(i,j) + d->depth(i+1,j));
    
    Ds(i,j) = MAX(0.00005, Ds(i,j));
    Dn(i,j) = MAX(0.00005, Dn(i,j));
    }
    
    SLICELOOP2
    {
    De(i,j) = ETAe(i,j)  + 0.5*(d->depth(i,j) + d->depth(i,j-1));
    Dw(i,j) = ETAw(i,j)  + 0.5*(d->depth(i,j) + d->depth(i,j+1));
    
    De(i,j) = MAX(0.00005, De(i,j));
    Dw(i,j) = MAX(0.00005, Dw(i,j));
    }
    
    // reconstruct U 
    precon->reconstruct_3D_x(p, pgc, d, UVEL, Us, Un);
    precon->reconstruct_3D_y(p, pgc, d, UVEL, Ue, Uw);
    precon->reconstruct_3D_z(p, pgc, d, UVEL, Ub, Ut);
    
    // reconstruct  V
    precon->reconstruct_3D_x(p, pgc, d, VVEL, Vs, Vn);
    precon->reconstruct_3D_y(p, pgc, d, VVEL, Ve, Vw);
    precon->reconstruct_3D_z(p, pgc, d, VVEL, Vb, Vt);
    
    // reconstruct  W
    precon->reconstruct_3D_x(p, pgc, d, WVEL, Ws, Wn);
    precon->reconstruct_3D_y(p, pgc, d, WVEL, We, Ww);
    precon->reconstruct_3D_z(p, pgc, d, WVEL, Wb, Wt);
    
    // signal speed
    pss->signal_speed_update(p, pgc, d, Us, Un, Ve, Vw, Ds, Dn, De, Dw);
}

void nhflow_HLL::start(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, slice &eta)
{
        if(ipol==1)
        aij_U(p,d,F,1,U,V,W);

        if(ipol==2 && p->j_dir==1)
        aij_V(p,d,F,2,U,V,W);

        if(ipol==3)
        aij_W(p,d,F,3,U,V,W);
}

double nhflow_HLL::aij_U(lexer* p,fdm_nhf* d, double *F, int ipol, double *UVEL, double *VVEL, double *WVEL)
{
    // flux x-dir
    ULOOP
    {
    Fs[IJK] = Us[IJK]*Us[IJK] 
            + (1.0/(Ds(i,j)))*(0.5*fabs(p->W22)*ETAs(i,j)*ETAs(i,j) + fabs(p->W22)*ETAs(i,j)*0.5*(d->depth(i,j) + d->depth(i-1,j)));
    
    Fn[IJK] = Un[IJK]*Un[IJK] 
            + (1.0/(Dn(i,j)))*(0.5*fabs(p->W22)*ETAn(i,j)*ETAn(i,j) + fabs(p->W22)*ETAn(i,j)*0.5*(d->depth(i,j) + d->depth(i+1,j)));
    }
    
    // flux y-dir
    VLOOP
    {
    Fe[IJK] = De(i,j)*Ue[IJK]*Ve[IJK];
    
    Fw[IJK] = Dw(i,j)*Uw[IJK]*Vw[IJK];
    }
    
    // flux z-dir
    WLOOP
    Fz[IJK] = 0.5*(d->omegaF[FIJK]*(Ub[IJK] + Ut[IJK]) - fabs(d->omegaF[FIJK])*(Ub[IJK] - Ut[IJK]));


    // HLL flux x-dir
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fx[IJK] = (d->Sn[IJK]*Fs[IJK] - d->Ss[IJK]*Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(Dn(i,j) - Ds(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fy[IJK] = (d->Sw[IJK]*Fe[IJK] - d->Se[IJK]*Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(Dw(i,j) - De(i,j)))/denom;
        }
    }
    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,11);
    pgc->start3V(p,Fz,12);
    
    LOOP
    {
    d->F[IJK] -= ((Fx[IJK] - Fx[Im1JK])/p->DXN[IP] 
                             + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir)
                             + (Fz[IJK] - Fz[IJKm1])/p->DZN[KP];
    }    
}

double nhflow_HLL::aij_V(lexer* p,fdm_nhf* d, double *F, int ipol, double *UVEL, double *VVEL, double *WVEL)
{
    // flux x-dir
    ULOOP
    {
    Fs[IJK] = Ds(i,j)*Us[IJK]*Vs[IJK];
    
    Fn[IJK] = Dn(i,j)*Un[IJK]*Vn[IJK];
    }
    
    // flux y-dir
    VLOOP
    {
    Fe[IJK] = De(i,j)*Ve[IJK]*Ve[IJK] 
            + 0.5*fabs(p->W22)*ETAe(i,j)*ETAe(i,j) + fabs(p->W22)*ETAe(i,j)*0.5*(d->depth(i,j) + d->depth(i,j-1));
    
    Fw[IJK] = Dw(i,j)*Vw[IJK]*Vw[IJK] 
            + 0.5*fabs(p->W22)*ETAw(i,j)*ETAw(i,j) + fabs(p->W22)*ETAw(i,j)*0.5*(d->depth(i,j) + d->depth(i,j+1));
    }
    
    // flux z-dir
    WLOOP
    Fz[IJK] = 0.5*(d->omegaF[FIJK]*(Vb[IJK] + Vt[IJK]) - fabs(d->omegaF[FIJK])*(Vb[IJK] - Vt[IJK]));


    // HLL flux x-dir
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fx[IJK] = (d->Sn[IJK]*Fs[IJK] - d->Ss[IJK]*Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(Dn(i,j) - Ds(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fy[IJK] = (d->Sw[IJK]*Fe[IJK] - d->Se[IJK]*Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(Dw(i,j) - De(i,j)))/denom;
        }
    }
    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,11);
    pgc->start3V(p,Fz,12);
    
    LOOP
    {
    d->G[IJK] -= (1.0/WLVL)*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP] 
                           + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir)
                           + (Fz[IJK] - Fz[IJKm1])/p->DZN[KP];
    }    
}

double nhflow_HLL::aij_W(lexer* p,fdm_nhf* d, double *F, int ipol, double *UVEL, double *VVEL, double *WVEL)
{
    // flux x-dir
    ULOOP
    {
    Fs[IJK] = Ds(i,j)*Us[IJK]*Ws[IJK];
    
    Fn[IJK] = Dn(i,j)*Un[IJK]*Wn[IJK];
    }
    
    // flux y-dir
    VLOOP
    {
    Fe[IJK] = De(i,j)*Ve[IJK]*We[IJK];
    
    Fw[IJK] = Dw(i,j)*Vw[IJK]*Ww[IJK];
    }
    
    // flux z-dir
    WLOOP
    Fz[IJK] = 0.5*(d->omegaF[FIJK]*(Wb[IJK] + Wt[IJK]) - fabs(d->omegaF[FIJK])*(Wb[IJK] - Wt[IJK]));


    // HLL flux x-dir
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fx[IJK] = (d->Sn[IJK]*Fs[IJK] - d->Ss[IJK]*Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(Dn(i,j) - Ds(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fy[IJK] = (d->Sw[IJK]*Fe[IJK] - d->Se[IJK]*Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(Dw(i,j) - De(i,j)))/denom;
        }
    }
    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,11);
    pgc->start3V(p,Fz,12);
    
    LOOP
    {
    d->H[IJK] -= (1.0/WLVL)*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP] 
                           + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir)
                           + (Fz[IJK] - Fz[IJKm1])/p->DZN[KP];
    }    
}



double nhflow_HLL::aij(lexer* p,fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, double *DX,double *DY, double *DZ)
{
    return 0.0;
}

