/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"idiff2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void idiff2::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, double sig, double alpha)
{
    count=0;

    sqd = (1.0/(p->dx*p->dx));

	LOOP
	{
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=visc(i,j,k);

	
//   M
	
	a->M.p[count]  +=   0.5*sqd*(vft*visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + vft*visc_ijk+ev_ijk/sig)
					+   0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)
					+   0.5*sqd*(vft*visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + vft*visc_ijk+ev_ijk/sig)
					+   0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)
					+   0.5*sqd*(vft*visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + vft*visc_ijk+ev_ijk/sig)
					+   0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig);

	 
	 a->M.s[count] -= 0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig);
	 a->M.n[count] -= 0.5*sqd*(vft*visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + vft*visc_ijk+ev_ijk/sig);
	 
	 a->M.e[count] -= 0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig);
	 a->M.w[count] -= 0.5*sqd*(vft*visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + vft*visc_ijk+ev_ijk/sig);
	 
	 a->M.b[count] -= 0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig);
	 a->M.t[count] -= 0.5*sqd*(vft*visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + vft*visc_ijk+ev_ijk/sig);
	 
	 ++count;
	}
}

void idiff2::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, double sig, double alpha)
{
}

