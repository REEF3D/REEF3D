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

#include"idiff2_FS_v2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

void idiff2_FS_v2::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc;
    
    pgc->start1(p,u,gcval_udiff);
	pgc->start2(p,v,gcval_vdiff);
	pgc->start3(p,w,gcval_wdiff);
    
	 
	count=0;
    if(p->j_dir==1)
    {
    VLOOP
    {
	visc = 0.5*(a->visc(i,j,k) + a->visc(i,j+1,k)) + 0.5*(a->eddyv(i,j,k) + a->eddyv(i,j+1,k));
    
	a->M.p[count] = 2.0*visc/(p->DYN[JP1]*p->DYP[JP])
				  + 2.0*visc/(p->DYN[JP]*p->DYP[JP])
				  + visc/(p->DXP[IP]*p->DXN[IP])
				  + visc/(p->DXP[IM1]*p->DXN[IP])
				  + visc/(p->DZP[KP]*p->DZN[KP])
				  + visc/(p->DZP[KM1]*p->DZN[KP])
				  + CPOR2/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((u(i,j+1,k)-u(i,j,k))*visc - (u(i-1,j+1,k)-u(i-1,j,k))*visc)/(p->DYP[JP]*p->DXN[IP])
						+  ((w(i,j+1,k)-w(i,j,k))*visc - (w(i,j+1,k-1)-w(i,j,k-1))*visc)/(p->DYP[JP]*p->DZN[KP])
									
						+ (CPOR2*v(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] = -visc/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -2.0*visc/(p->DYN[JP]*p->DYP[JP]);
	 a->M.w[count] = -2.0*visc/(p->DYN[JP1]*p->DYP[JP]);
	 
	 a->M.b[count] = -visc/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc/(p->DZP[KP]*p->DZN[KP]);
	 
	 ++count;
	}
    
    n=0;
	VLOOP
	{
		if(p->flag2[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*v(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag2[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*v(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag2[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*v(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag2[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*v(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag2[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*v(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag2[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*v(i,j,k+1);
		a->M.t[n] = 0.0;
		}

	++n;
	}


	psolv->start(p,a,pgc,v,a->rhsvec,2);
    }

    pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);

	time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"vdiffiter: "<<p->viter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}





