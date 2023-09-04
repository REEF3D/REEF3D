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

#include"idiff2_FS.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

void idiff2_FS::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;

	count=0;
    if(p->k_dir==1)
    {
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
	
	visc_ddx_p = 0.25*(visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_ip_j_k+ev_ip_j_k + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1));
    
	visc_ddx_m = 0.25*(visc_im_j_k+ev_im_j_k + a->visc(i-1,j,k+1)+a->eddyv(i-1,j,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp);
    
	visc_ddy_p = 0.25*(visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_i_jp_k+ev_i_jp_k + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1));
    
	visc_ddy_m = 0.25*(visc_i_jm_k+ev_i_jm_k + a->visc(i,j-1,k+1)+a->eddyv(i,j-1,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp);
    
	a->M.p[count] = 2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				  + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				  + CPOR3/(alpha*p->dt);
				  
	a->rhsvec.V[count] +=  ((a->u(i,j,k+1)-u(i,j,k))*visc_ddx_p - (u(i-1,j,k+1)-u(i-1,j,k))*visc_ddx_m)/(p->DZP[KP]*p->DXN[IP])
						+  ((a->v(i,j,k+1)-v(i,j,k))*visc_ddy_p - (v(i,j-1,k+1)-v(i,j-1,k))*visc_ddy_m)/(p->DZP[KP]*p->DYN[JP])
									
						+ (CPOR3*w(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] = -visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] = -2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP]);
	 a->M.t[count] = -2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP]);
    
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
    
	pgc->start3(p,w,gcval_w);
	
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}


void idiff2_FS::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &w_in, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;

	count=0;
    
    WLOOP
    diff(i,j,k) = w_in(i,j,k);
    
	pgc->start3(p,diff,gcval_w);

	count=0;
    if(p->k_dir==1)
    {
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
	
	visc_ddx_p = 0.25*(visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_ip_j_k+ev_ip_j_k + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1));
    
	visc_ddx_m = 0.25*(visc_im_j_k+ev_im_j_k + a->visc(i-1,j,k+1)+a->eddyv(i-1,j,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp);
    
	visc_ddy_p = 0.25*(visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_i_jp_k+ev_i_jp_k + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1));
    
	visc_ddy_m = 0.25*(visc_i_jm_k+ev_i_jm_k + a->visc(i,j-1,k+1)+a->eddyv(i,j-1,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp);
    
	a->M.p[count] = 2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				  + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				  + CPOR3/(alpha*p->dt);
				  
	a->rhsvec.V[count] +=  ((a->u(i,j,k+1)-u(i,j,k))*visc_ddx_p - (u(i-1,j,k+1)-u(i-1,j,k))*visc_ddx_m)/(p->DZP[KP]*p->DXN[IP])
						+  ((a->v(i,j,k+1)-v(i,j,k))*visc_ddy_p - (v(i,j-1,k+1)-v(i,j-1,k))*visc_ddy_m)/(p->DZP[KP]*p->DYN[JP])
									
						+ (CPOR3*w_in(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] = -visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] = -2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP]);
	 a->M.t[count] = -2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP]);
    
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
    
	
	psolv->start(p,a,pgc,diff,a->rhsvec,3);
    }
    
    
	pgc->start3(p,diff,gcval_w);
	
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}

