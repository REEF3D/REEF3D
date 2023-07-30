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

void idiff2_FS::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;

	count=0;
    if(p->j_dir==1)
    {
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
	
	visc_ddx_p = 0.25*(visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k + visc_ip_j_k+ev_ip_j_k + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k));
    
	visc_ddx_m = 0.25*(visc_im_j_k+ev_im_j_k + a->visc(i-1,j+1,k)+a->eddyv(i-1,j+1,k) + visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k);
    
	visc_ddz_p = 0.25*(visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k + visc_i_j_kp+ev_i_j_kp + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1));

	visc_ddz_m = 0.25*(visc_i_j_km+ev_i_j_km + a->visc(i,j+1,k-1)+a->eddyv(i,j+1,k-1) + visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k);
    
	
	a->M.p[count] = 2.0*(visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP1]*p->DYP[JP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DYN[JP]*p->DYP[JP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				  + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				  + CPOR2/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((u(i,j+1,k)-u(i,j,k))*visc_ddx_p - (u(i-1,j+1,k)-u(i-1,j,k))*visc_ddx_m)/(p->DYP[JP]*p->DXN[IP])
						+  ((w(i,j+1,k)-w(i,j,k))*visc_ddz_p - (w(i,j+1,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DYP[JP]*p->DZN[KP])
									
						+ (CPOR2*v(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -2.0*(visc_ijk+ev_ijk)/(p->DYN[JP]*p->DYP[JP]);
	 a->M.w[count] = -2.0*(visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP1]*p->DYP[JP]);
	 
	 a->M.b[count] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
     
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
    
	pgc->start2(p,v,gcval_v);

    time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"vdiffiter: "<<p->viter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}


void idiff2_FS::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &v_in, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;
    
    VLOOP
    diff(i,j,k) = v_in(i,j,k);
    
    pgc->start2(p,diff,gcval_v);

	count=0;

	count=0;
    if(p->j_dir==1)
    {
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
	
	visc_ddx_p = 0.25*(visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k + visc_ip_j_k+ev_ip_j_k + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k));
    
	visc_ddx_m = 0.25*(visc_im_j_k+ev_im_j_k + a->visc(i-1,j+1,k)+a->eddyv(i-1,j+1,k) + visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k);
    
	visc_ddz_p = 0.25*(visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k + visc_i_j_kp+ev_i_j_kp + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1));

	visc_ddz_m = 0.25*(visc_i_j_km+ev_i_j_km + a->visc(i,j+1,k-1)+a->eddyv(i,j+1,k-1) + visc_ijk+ev_ijk + visc_i_jp_k+ev_i_jp_k);
    
	
	a->M.p[count] = 2.0*(visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP1]*p->DYP[JP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DYN[JP]*p->DYP[JP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				  + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				  + CPOR2/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((u(i,j+1,k)-u(i,j,k))*visc_ddx_p - (u(i-1,j+1,k)-u(i-1,j,k))*visc_ddx_m)/(p->DYP[JP]*p->DXN[IP])
						+  ((w(i,j+1,k)-w(i,j,k))*visc_ddz_p - (w(i,j+1,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DYP[JP]*p->DZN[KP])
									
						+ (CPOR2*v_in(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[count] = -2.0*(visc_ijk+ev_ijk)/(p->DYN[JP]*p->DYP[JP]);
	 a->M.w[count] = -2.0*(visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP1]*p->DYP[JP]);
	 
	 a->M.b[count] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
     
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


	psolv->start(p,a,pgc,diff,a->rhsvec,2);
    }
    
	pgc->start2(p,diff,gcval_v);
	
    time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"vdiffiter: "<<p->viter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}



