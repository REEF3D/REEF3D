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

#include"idiff2_FS_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

void idiff2_FS_2D::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field& b, field& visc, double sig, double alpha)
{
    starttime=pgc->timer();
    
    sqd = (1.0/(p->dx*p->dx));
    
    
    count=0;
	LOOP
	{
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=visc(i,j,k);

	
//   M
	
	a->M.p[count]  =    0.5*sqd*(vft*visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + vft*visc_ijk+ev_ijk/sig)
					+   0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig)
					+   0.5*sqd*(vft*visc(i,j+1,k)+a->eddyv(i,j+1,k)/sig + vft*visc_ijk+ev_ijk/sig)
					+   0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j-1,k)+a->eddyv(i,j-1,k)/sig)
					+   0.5*sqd*(vft*visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + vft*visc_ijk+ev_ijk/sig)
					+   0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig)
					+   1.0/(alpha*p->dt);
    
    a->rhsvec.V[count] += b(i,j,k)/(alpha*p->dt);

	 
	 a->M.s[count] = -0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i-1,j,k)+a->eddyv(i-1,j,k)/sig);
	 a->M.n[count] = -0.5*sqd*(vft*visc(i+1,j,k)+a->eddyv(i+1,j,k)/sig + vft*visc_ijk+ev_ijk/sig);
	 
	 a->M.b[count] = -0.5*sqd*(vft*visc_ijk+ev_ijk/sig + vft*visc(i,j,k-1)+a->eddyv(i,j,k-1)/sig);
	 a->M.t[count] = -0.5*sqd*(vft*visc(i,j,k+1)+a->eddyv(i,j,k+1)/sig + vft*visc_ijk+ev_ijk/sig);
	 
	 ++count;
	}
    
    psolv->start(p,a,pgc,b,a->xvec,a->rhsvec,4,1,p->D29);
    time=pgc->timer()-starttime;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"scalar_diffiter: "<<p->solveriter<<"  scalar_difftime: "<<setprecision(3)<<time<<endl;
}

void idiff2_FS_2D::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field& b, field& visc, double sig, double alpha)
{
}
