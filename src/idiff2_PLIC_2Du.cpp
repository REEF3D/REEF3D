/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"idiff2_PLIC_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

idiff2_PLIC_2D::idiff2_PLIC_2D(lexer* p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

idiff2_PLIC_2D::~idiff2_PLIC_2D()
{
}

void idiff2_PLIC_2D::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &u_in, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	double visc_ddz_p,visc_ddz_m, H_ddz_p, H_ddz_m;
    
    ULOOP
    diff(i,j,k) = u_in(i,j,k);
    
    pgc->start1(p,diff,gcval_u);


	count=0;
    if(p->i_dir==1)
    {
	ULOOP
	{
    
    if(p->F92==22 || p->F92==32)
    {
        ev_ijk=a->vof(i,j,k)*a->eddyv(i,j,k);
        ev_ip_j_k=a->vof(i+1,j,k)*a->eddyv(i+1,j,k);
        ev_i_j_km=a->vof(i,j,k-1)*a->eddyv(i,j,k-1);
        ev_i_j_kp=a->vof(i,j,k+1)*a->eddyv(i,j,k+1);
        ev_ip_j_kp=a->vof(i+1,j,k+1)*a->eddyv(i+1,j,k+1);
        ev_ip_j_km=a->vof(i+1,j,k-1)*a->eddyv(i+1,j,k-1);
    }
    else
    {
        ev_ijk=a->eddyv(i,j,k);
        ev_ip_j_k=a->eddyv(i+1,j,k);
        ev_i_j_km=a->eddyv(i,j,k-1);
        ev_i_j_kp=a->eddyv(i,j,k+1);
        ev_ip_j_kp=a->eddyv(i+1,j,k+1);
        ev_ip_j_km=a->eddyv(i+1,j,k-1);
    }
	visc_ijk=a->visc(i,j,k);
	visc_ip_j_k=a->visc(i+1,j,k);
	visc_i_j_km=a->visc(i,j,k-1);
	visc_i_j_kp=a->visc(i,j,k+1);
    
    if(p->F92==3||p->F92==32)
    {
        H_ddz_p=(a->vof_nt(i,j,k)+a->vof_st(i+1,j,k)+a->vof_nb(i,j,k+1)+a->vof_sb(i+1,j,k+1))*0.25;
        H_ddz_m=(a->vof_nb(i,j,k)+a->vof_sb(i+1,j,k)+a->vof_nt(i,j,k-1)+a->vof_st(i+1,j,k-1))*0.25;
        
        visc_ddz_p = H_ddz_p*p->W2+(1.0-H_ddz_p)*p->W4+(ev_ijk+ev_ip_j_k+ev_i_j_kp+ev_ip_j_kp)*0.25;
        visc_ddz_m = H_ddz_m*p->W2+(1.0-H_ddz_m)*p->W4+(ev_ijk+ev_ip_j_k+ev_i_j_km+ev_ip_j_km)*0.25;
	}
    else
    {
        visc_ddz_p = (visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k + visc_i_j_kp+ev_i_j_kp + a->visc(i+1,j,k+1)+ev_ip_j_kp)*0.25;
        visc_ddz_m = (a->visc(i,j,k-1)+ev_i_j_km + a->visc(i+1,j,k-1)+ev_ip_j_km + visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k)*0.25;
    }
        
	a->M.p[count] =  2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP])
				   + 2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP])
				   + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				   + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				   + CPOR1/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((w(i+1,j,k)-w(i,j,k))*visc_ddz_p - (w(i+1,j,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DXP[IP]*p->DZN[KP])

						 + (CPOR1*u_in(i,j,k))/(alpha*p->dt);

	 
	 a->M.s[count] = -2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP]);
	 a->M.n[count] = -2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP]);
	 
	 a->M.b[count] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
	 
	 ++count;
	}
    
    n=0;
	ULOOP
	{
		if(p->flag1[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*u(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag1[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*u(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag1[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*u(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag1[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*u(i,j,k+1);
		a->M.t[n] = 0.0;
		}

	++n;
	}
	
	psolv->start(p,a,pgc,diff,a->rhsvec,1);
    }
	
    pgc->start1(p,diff,gcval_u);

    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
}


