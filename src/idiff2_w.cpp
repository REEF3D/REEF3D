/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
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

void idiff2::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
     count=0;
	 
	 double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;

	WLOOP
    {
	ev_ijk=a->eddyv(i,j,k);
	ev_im_j_k=a->eddyv(i-1,j,k);
	ev_ip_j_k=a->eddyv(i+1,j,k);
	ev_i_jm_k=a->eddyv(i,j-1,k);
	ev_i_jp_k=a->eddyv(i,j+1,k);
	ev_i_j_kp=a->eddyv(i,j,k+1);
	
	visc_ijk=a->visc(i,j,k);
	visc_im_j_k=a->visc(i-1,j,k);
	visc_ip_j_k=a->visc(i+1,j,k);
	visc_i_jm_k=a->visc(i,j-1,k);
	visc_i_jp_k=a->visc(i,j+1,k);
	visc_i_j_kp=a->visc(i,j,k+1);
	
	visc_ddx_p = (visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_ip_j_k+ev_ip_j_k + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddx_m = (visc_im_j_k+ev_im_j_k + a->visc(i-1,j,k+1)+a->eddyv(i-1,j,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp)*0.25;
	visc_ddy_p = (visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_i_jp_k+ev_i_jp_k + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
	visc_ddy_m = (visc_i_jm_k+ev_i_jm_k + a->visc(i,j-1,k+1)+a->eddyv(i,j-1,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp)*0.25;
	
    
	a->M.p[count] += 2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				  + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				  + CPOR3/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((u(i,j,k+1)-u(i,j,k))*visc_ddx_p - (u(i-1,j,k+1)-u(i-1,j,k))*visc_ddx_m)/(p->DZP[KP]*p->DXN[IP])


						+  ((v(i,j,k+1)-v(i,j,k))*visc_ddy_p - (v(i,j-1,k+1)-v(i,j-1,k))*visc_ddy_m)/(p->DZP[KP]*p->DYN[JP])
									
						+ (CPOR3*w(i,j,k))/(alpha*p->dt);
									
	 a->M.s[count] -= visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] -= visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] -= visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] -= visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] -= 2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP]);
	 a->M.t[count] -= 2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP]);
	 
	 ++count;
	}

}


