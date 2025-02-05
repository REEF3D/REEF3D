/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_idiff_2D.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"

void nhflow_idiff_2D::diff_scalar(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *F, double sig, double alpha)
{
	n=0;

	LOOP
	{
	visc = d->VISC[IJK] + d->EV[IJK]/sig;
    
    sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
	
//   M
	
	d->M.p[n]  +=        visc/(p->DXN[IP]*p->DXP[IM1])
					+   visc/(p->DXN[IP]*p->DXP[IP])
					+   visc/(p->DYN[JP]*p->DYP[JM1])*p->y_dir
					+   visc/(p->DYN[JP]*p->DYP[JP])*p->y_dir
					+   (visc*sigxyz2)/(p->DZN[KP]*p->DZP[KM1])
					+   (visc*sigxyz2)/(p->DZN[KP]*p->DZP[KP]);
    
	 d->M.s[n] -= visc/(p->DXN[IP]*p->DXP[IM1]);
	 d->M.n[n] -= visc/(p->DXN[IP]*p->DXP[IP]);
	 
	 d->M.b[n] -= (visc*sigxyz2)/(p->DZN[KP]*p->DZP[KM1]);
	 d->M.t[n] -= (visc*sigxyz2)/(p->DZN[KP]*p->DZP[KP]);
     
     d->rhsvec.V[n] +=       visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(F[Ip1JKp1] - F[Im1JKp1] - F[Ip1JKm1] + F[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]));
	 
	 ++n;
	}
}
