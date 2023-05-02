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

#include"nhflow_flux_FOU.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_flux_FOU::nhflow_flux_FOU(lexer* p, patchBC_interface *ppBC) : nhflow_flux_reconstruct(p,ppBC)
                                                                    
{
    pBC = ppBC;
    
    p->Darray(Fs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fw,p->imax*p->jmax*(p->kmax+2));
}

nhflow_flux_FOU::~nhflow_flux_FOU()
{
}

void nhflow_flux_FOU::face_flux_2D(lexer* p, fdm_nhf*, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
{
}

void nhflow_flux_FOU::face_flux_3D(lexer *p, ghostcell *pgc, fdm_nhf *d, slice & eta, double *U, double *V, double *Fx, double *Fy)
{
    double Ss,Sn,Se,Sw;
    double USx,USy;
    double Ds,Dn,De,Dw;
    double DSx,DSy;
    double denom;
    
    // FOU flux
    ULOOP
    {
    Ss = 0.5*(U[IJK] + U[Ip1JK]);

        // final flux x-dir
        if(Ss>=0.0)
        Fx[IJK] = U[IJK];
        
        if(Ss<0.0)
        Fx[IJK] = U[Ip1JK];
    }
    
    
    VLOOP
    {
    Ss = 0.5*(V[IJK] + V[IJp1K]);

        // final flux y-dir
        if(Ss>=0.0)
        Fy[IJK] = V[IJK];
        
        if(Ss<0.0)
        Fy[IJK] = V[IJp1K];
    }
    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,10);
	
}


