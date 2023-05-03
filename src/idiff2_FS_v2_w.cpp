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

void idiff2_FS_v2::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc;
    
	count=0;
    
    pgc->start1(p,u,gcval_udiff);
	pgc->start2(p,v,gcval_vdiff);
	pgc->start3(p,w,gcval_wdiff);

	count=0;
    if(p->k_dir==1)
    {
    WLOOP
    {
    visc = 0.5*(a->visc(i,j,k) + a->visc(i,j,k+1)) + 0.5*(a->eddyv(i,j,k) + a->eddyv(i,j,k+1));
    
	a->M.p[count] = 2.0*visc/(p->DZN[KP1]*p->DZP[KP])
				  + 2.0*visc/(p->DZN[KP]*p->DZP[KP])
				  + visc/(p->DXP[IP]*p->DXN[IP])
				  + visc/(p->DXP[IM1]*p->DXN[IP])
				  + visc/(p->DYP[JP]*p->DYN[JP])
				  + visc/(p->DYP[JM1]*p->DYN[JP])
				  + CPOR3/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((a->u(i,j,k+1)-u(i,j,k))*visc - (u(i-1,j,k+1)-u(i-1,j,k))*visc)/(p->DZP[KP]*p->DXN[IP])
						+  ((a->v(i,j,k+1)-v(i,j,k))*visc - (v(i,j-1,k+1)-v(i,j-1,k))*visc)/(p->DZP[KP]*p->DYN[JP])
									
						+ (CPOR3*w(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] = -visc/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -visc/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] = -visc/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] = -2.0*visc/(p->DZN[KP]*p->DZP[KP]);
	 a->M.t[count] = -2.0*visc/(p->DZN[KP1]*p->DZP[KP]);
	 
	 ++count;
	}
    
    n=0;
    WLOOP
	{
		if(p->flag3[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*w(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag3[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*w(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag3[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*w(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag3[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*w(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag3[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*w(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag3[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*w(i,j,k+1);
		a->M.t[n] = 0.0;
		}

	++n;
	}
    
	
	psolv->start(p,a,pgc,w,a->rhsvec,3);
    }
    
    
    pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}



