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
    double uval, vval, wval;
    double Cf=0.0001;
    
    if(p->X38==1)
    {
    xc = xp + p->X43*nx*p->DXP[IP];
    yc = yp + p->X43*ny*p->DYP[JP];
    zc = zp + p->X43*nz*p->DZP[KP]*WL(i,j);
        
    uval   = p->ccipol4V(d->U, WL, d->bed, xc, yc, zc); 
    vval   = p->ccipol4V(d->V, WL, d->bed, xc, yc, zc); 
    wval   = p->ccipol4V(d->W, WL, d->bed, xc, yc, zc); 
    

	Fv_x = 0.5*p->W1*Cf*uval*fabs(uval);
    Fv_y = 0.5*p->W1*Cf*vval*fabs(vval);
    Fv_z = 0.5*p->W1*Cf*wval*fabs(wval);
    
    if(p->j_dir==0)
    Fv_y = 0.0;
    }
    
    
}