/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"6DOF_obj.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::hydrodynamic_viscous_forces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL,
                                            double &Fv_x, double &Fv_y, double &Fv_z, double A_triang,
                                            double xp, double yp, double zp, double nx, double ny, double nz)
{
    double uval, vval, wval, Uabs;
    double ReL, CF0, Cf;

    Fv_x = Fv_y = Fv_z = 0.0;

    // ITTC-1957 friction line combined with a form factor (1+k):
    // R_F = 0.5*rho*(1+k)*CF0*U^2*S , with CF0 = 0.075/(log10(Re_L)-2)^2
    // applied locally as a quadratic wall shear stress using the flow
    // velocity sampled next to the hull, oriented along the local
    // relative velocity direction (X39_Lwl: reference/waterline length,
    // X39_k: form factor k).
    if(p->X39==1)
    {
    xc = xp + p->X43*nx*p->DXP[IP];
    yc = yp + p->X43*ny*p->DYP[JP];
    zc = zp + p->X43*nz*p->DZP[KP]*WL(i,j);

    uval   = p->ccipol4V(d->U, WL, d->bed, xc, yc, zc);
    vval   = p->ccipol4V(d->V, WL, d->bed, xc, yc, zc);
    wval   = p->ccipol4V(d->W, WL, d->bed, xc, yc, zc);

    Uabs = sqrt(uval*uval + vval*vval + wval*wval);

    ReL = Uabs*p->X39_Lwl/p->W2;

    CF0 = 0.0;
    if(ReL>1.0e4)
    CF0 = 0.075/pow(log10(ReL)-2.0,2.0);

    Cf = (1.0 + p->X39_k)*CF0;

	Fv_x = 0.5*p->W1*Cf*uval*fabs(uval)*A_triang;
    Fv_y = 0.5*p->W1*Cf*vval*fabs(vval)*A_triang;
    Fv_z = 0.5*p->W1*Cf*wval*fabs(wval)*A_triang;

    if(p->j_dir==0)
    Fv_y = 0.0;
    }


}