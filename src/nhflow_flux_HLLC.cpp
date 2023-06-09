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

#include"nhflow_flux_HLLC.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_reconstruct_WENO.h"

nhflow_flux_HLLC::nhflow_flux_HLLC(lexer* p, patchBC_interface *ppBC) : ETAs(p),ETAn(p),ETAe(p),ETAw(p)
{
    pBC = ppBC;
    
    if(p->A543==2)
    precon = new nhflow_reconstruct_hires(p,ppBC);
    
    if(p->A543==4)
    precon = new nhflow_reconstruct_weno(p,ppBC);
    
    p->Darray(DU,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DV,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Fs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fw,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Us,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Un,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ve,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vw,p->imax*p->jmax*(p->kmax+2));
}

nhflow_flux_HLLC::~nhflow_flux_HLLC()
{
}

void nhflow_flux_HLLC::face_flux_2D(lexer* p, fdm_nhf*, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
{
}

void nhflow_flux_HLLC::face_flux_3D(lexer *p, ghostcell *pgc, fdm_nhf *d, slice & eta, double *U, double *V, double *Fx, double *Fy)
{
    double Ss,Sn,Se,Sw,SS;
    double USx,USy;
    double Ds,Dn,De,Dw;
    double DSx,DSy;
    double denom;
    double FsS,FnS;
    double FeS,FwS;
    
    // velocity depth
    LOOP
    {
    DU[IJK] = U[IJK]*d->WL[IJ];
    DV[IJK] = V[IJK]*d->WL[IJ];
    }
    
    pgc->start1V(p,DU,10);
    pgc->start2V(p,DV,11);
    
    // reconstruct eta
    precon->reconstruct_2D(p, pgc, d, eta, ETAs, ETAn, ETAe, ETAw);
    
    // reconstruct U and V
    precon->reconstruct_3D(p, pgc, d, DU, DV, Fs, Fn, Fe, Fw);
    precon->reconstruct_3D(p, pgc, d, U, V, Us, Un, Ve, Vw);
    
    // HLLC flux
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
    SS = USx;
    
    if(p->wet[Ip1J]==0)
    {
    Ss = Us[IJK] - sqrt(9.81*Ds);
    Sn = Us[IJK] + 2.0*sqrt(9.81*Ds);
    SS = Ss;
    }
    
    if(p->wet[Im1J]==0)
    {
    Ss = Un[IJK] - 2.0*sqrt(9.81*Dn);
    Sn = Un[IJK] + sqrt(9.81*Dn);
    SS = Sn;
    }
    
    if(p->wet[IJ]==0)
    {
    Ss=Sn=0.0;
    USx=0.0;
    }
    
    FsS = Ds*(Ss - Us[IJK] + 1.0e-10)/(Ss - SS + 1.0e-10);
    FnS = Dn*(Sn - Un[IJK] + 1.0e-10)/(Sn - SS + 1.0e-10);
    
        // final flux x-dir
        if(Ss>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(Sn<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        if(SS>=0.0)
        Fx[IJK] = Fs[IJK] + Ss*(FsS - Us[IJK]);
        
        else
        Fx[IJK] = Fn[IJK] + Sn*(FnS - Un[IJK]);
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
    SS = USy;
    
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


    FeS = De*(Se - Ve[IJK] + 1.0e-10)/(Se - SS + 1.0e-10);
    FwS = Dw*(Sw - Vw[IJK] + 1.0e-10)/(Sw - SS + 1.0e-10);
 
        // final flux x-dir
        if(Ss>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(Sn<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        if(SS>=0.0)
        Fy[IJK] = Fe[IJK] + Ss*(FeS - Ve[IJK]);
        
        else
        Fy[IJK] = Fw[IJK] + Sn*(FwS - Vw[IJK]);
    }
    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,10);
	
}


