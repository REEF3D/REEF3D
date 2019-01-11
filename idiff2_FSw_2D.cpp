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

void idiff2_FS_2D::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;
    
    if(p->D24==0)
    alpha=1.0;
	
	count=0;
    
    pgc->start3(p,w,gcval_w);

    sqd = (1.0/(p->dx*p->dx));


	count=0;

    WLOOP
    {
	ev_ijk=a->eddyv(i,j,k);
	ev_im_j_k=a->eddyv(i-1,j,k);
	ev_ip_j_k=a->eddyv(i+1,j,k);
	ev_i_j_kp=a->eddyv(i,j,k+1);
	
	visc_ijk=a->visc(i,j,k);
	visc_im_j_k=a->visc(i-1,j,k);
	visc_ip_j_k=a->visc(i+1,j,k);
	visc_i_j_kp=a->visc(i,j,k+1);
	
	visc_ddx_p = (vfm*visc_ijk+ev_ijk + vfm*visc_i_j_kp+ev_i_j_kp + vfm*visc_ip_j_k+ev_ip_j_k + vfm*a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddx_m = (vfm*visc_im_j_k+ev_im_j_k + vfm*a->visc(i-1,j,k-1)+a->eddyv(i-1,j,k-1) + vfm*visc_ijk+ev_ijk + vfm*visc_i_j_kp+ev_i_j_kp)*0.25;
    
	a->M.p[count] = 2.0*(vfm*visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP]*p->DZP[KP])
				  + 2.0*(vfm*visc_ijk+ev_ijk)/(p->DZN[KM1]*p->DZP[KP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + CPOR3/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((u(i,j,k+1)-u(i,j,k))*visc_ddx_p - (u(i-1,j,k+1)-u(i-1,j,k))*visc_ddx_m)/(p->DZP[KP]*p->DXN[IP])
									
						+  a->M.p[count]*w(i,j,k)*(1.0/p->N54-1.0)
						+ (CPOR3*w(i,j,k))/(alpha*p->dt);
									
	a->M.p[count] /= p->N54;
	 
	 a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.b[count] = -2.0*(vfm*visc_ijk+ev_ijk)/(p->DZN[KM1]*p->DZP[KP]);
	 a->M.t[count] = -2.0*(vfm*visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP]*p->DZP[KP]);
	 
	 ++count;
	}
    
	
	psolv->start(p,a,pgc,w,a->xvec,a->rhsvec,3,gcval_w,p->D29);

    
    pgc->start3(p,w,gcval_w);
	
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}



