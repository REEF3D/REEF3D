/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"idiff2_FS.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

void idiff2_FS::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;
    
    if(p->D24==0)
    alpha=1.0;
    
    pgc->start2(p,v,gcval_v);
	 
	count=0;

    sqd = (1.0/(p->dx*p->dx));


	count=0;
    if(p->j_dir==1)
    VLOOP
    {
	ev_ijk=a->eddyv(i,j,k);
	ev_im_j_k=a->eddyv(i-1,j,k);
	ev_ip_j_k=a->eddyv(i+1,j,k);
	ev_i_jp_k=a->eddyv(i,j+1,k);
	ev_i_j_km=a->eddyv(i,j,k-1);
	ev_i_j_kp=a->eddyv(i,j,k+1);
	
	visc_ijk=a->visc(i,j,k);
	visc_im_j_k=a->visc(i-1,j,k);
	visc_ip_j_k=a->visc(i+1,j,k);
	visc_i_jp_k=a->visc(i,j+1,k);
	visc_i_j_km=a->visc(i,j,k-1);
	visc_i_j_kp=a->visc(i,j,k+1);
	
	visc_ddx_p = (vfm*visc_ijk+ev_ijk + vfm*visc_i_jp_k+ev_i_jp_k + vfm*visc_ip_j_k+ev_ip_j_k + vfm*a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
	visc_ddx_m = (vfm*visc_im_j_k+ev_im_j_k + vfm*a->visc(i-1,j+1,k)+a->eddyv(i-1,j+1,k) + vfm*visc_ijk+ev_ijk + vfm*visc_i_jp_k+ev_i_jp_k)*0.25;
	visc_ddz_p = (vfm*visc_ijk+ev_ijk + vfm*visc_i_jp_k+ev_i_jp_k + vfm*visc_i_j_kp+ev_i_j_kp + vfm*a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
	visc_ddz_m = (vfm*visc_i_j_km+ev_i_j_km + vfm*a->visc(i,j+1,k-1)+a->eddyv(i,j+1,k-1) + vfm*visc_ijk+ev_ijk + vfm*visc_i_jp_k+ev_i_jp_k)*0.25;
	
	
	a->M.p[count] = 2.0*(vfm*visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP]*p->DYP[JP])
				  + 2.0*(vfm*visc_ijk+ev_ijk)/(p->DYN[JM1]*p->DYP[JP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				  + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				  + CPOR2/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((u(i,j+1,k)-u(i,j,k))*visc_ddx_p - (u(i-1,j+1,k)-u(i-1,j,k))*visc_ddx_m)/(p->DYP[JP]*p->DXN[IP])

						+  ((w(i,j+1,k)-w(i,j,k))*visc_ddz_p - (w(i,j+1,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DYP[JP]*p->DZN[KP])
									
						+  a->M.p[count]*v(i,j,k)*(1.0/p->N54-1.0)
						+ (CPOR2*v(i,j,k))/(alpha*p->dt);
						
	a->M.p[count] /= p->N54;
	 
	 a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -2.0*(vfm*visc_ijk+ev_ijk)/(p->DYN[JM1]*p->DYP[JP]);
	 a->M.w[count] = -2.0*(vfm*visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP]*p->DYP[JP]);
	 
	 a->M.b[count] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
	 
	 ++count;
	}


	psolv->start(p,a,pgc,v,a->xvec,a->rhsvec,2,gcval_v,p->D29);
    
	time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
	cout<<"vdiffiter: "<<p->viter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}





