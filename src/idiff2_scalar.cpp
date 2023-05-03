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

#include"idiff2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void idiff2::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &eddyv, double sig, double alpha)
{
    count=0;

	LOOP
	{
	ev_ijk=eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);

	
//   M
	
	a->M.p[count]  +=    0.5*(a->visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + visc_ijk+ev_ijk/sig)/(p->DXN[IP]*p->DXP[IM1])
					+   0.5*(visc_ijk+ev_ijk/sig + a->visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)/(p->DXN[IP]*p->DXP[IP])
					+   0.5*(a->visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + visc_ijk+ev_ijk/sig)/(p->DYN[JP]*p->DYP[JM1])*p->y_dir
					+   0.5*(visc_ijk+ev_ijk/sig + a->visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)/(p->DYN[JP]*p->DYP[JP])*p->y_dir
					+   0.5*(a->visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + visc_ijk+ev_ijk/sig)/(p->DZN[KP]*p->DZP[KM1])
					+   0.5*(visc_ijk+ev_ijk/sig + a->visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)/(p->DZN[KP]*p->DZP[KP]);
    
	 a->M.s[count] -= 0.5*(visc_ijk+ev_ijk/sig + a->visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)/(p->DXN[IP]*p->DXP[IM1]);
	 a->M.n[count] -= 0.5*(a->visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + visc_ijk+ev_ijk/sig)/(p->DXN[IP]*p->DXP[IP]);
	 
	 a->M.e[count] -= 0.5*(visc_ijk+ev_ijk/sig + a->visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)/(p->DYN[JP]*p->DYP[JM1])*p->y_dir;
	 a->M.w[count] -= 0.5*(a->visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + visc_ijk+ev_ijk/sig)/(p->DYN[JP]*p->DYP[JP])*p->y_dir;
	 
	 a->M.b[count] -= 0.5*(visc_ijk+ev_ijk/sig + a->visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)/(p->DZN[KP]*p->DZP[KM1]);
	 a->M.t[count] -= 0.5*(a->visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + visc_ijk+ev_ijk/sig)/(p->DZN[KP]*p->DZP[KP]);
	 
	 ++count;
	}
}

void idiff2::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, field &eddyv, double sig, double alpha)
{
}

