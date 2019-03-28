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

#include"sflow_idiff.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_idiff::sflow_idiff(lexer* p)
{
}

sflow_idiff::~sflow_idiff()
{
}

void sflow_idiff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{
    double visc;
    
    
    
	SLICELOOP1
    {

	b->F(i,j) +=  (visc/(p->dx*p->dx))*
    
                (2.0*(u(i+1,j) - 2.0*u(i,j) + u(i-1,j))
                    +(u(i,j+1) - 2.0*u(i,j) + u(i,j-1))
                
                + (v(i+1,j)-v(i,j)) - (v(i+1,j-1)-v(i,j-1)));
                                        
    }
    
    SLICELOOP1
    {
	visc = p->W2 + 0.5*(b->eddyv(i,j) + b->eddyv(i+1,j));

        
	b->M.p[count] =  4.0*visc/(p->dx*p->dx);
                   
				   + 1.0/(alpha*p->dt);
				  
	b->rhsvec.V[count] += (visc/(p->dx*p->dx))*((v(i+1,j)-v(i,j)) - (v(i+1,j-1)-v(i,j-1)))
									
						 + b->M.p[count]*u(i,j)*(1.0/p->N54-1.0)
						 + (u(i,j))/(alpha*p->dt);
									
	 b->M.p[count] /= p->N54;
	 
	 b->M.s[count] = -2.0*visc/(p->dx*p->dx);
	 b->M.n[count] = -2.0*visc/(p->dx*p->dx);
	 
	 b->M.e[count] = -2.0*visc/(p->dx*p->dx);
	 b->M.w[count] = -2.0*visc/(p->dx*p->dx);
 
	 ++count;
	}
    
    /*
    starttime=pgc->timer();
	double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;
    
    if(p->D24==0)
    alpha=1.0;
    
    pgc->start1(p,u,gcval_u);

    count=0;

    sqd = (1.0/(p->dx*p->dx));
	
    
	count=0;
    if(p->i_dir==1)
    {
	ULOOP
	{
	ev_ijk=b->eddyv(i,j,k);
	ev_ip_j_k=b->eddyv(i+1,j,k);
	ev_i_jm_k=b->eddyv(i,j-1,k);
	ev_i_jp_k=b->eddyv(i,j+1,k);
	ev_i_j_km=b->eddyv(i,j,k-1);
	ev_i_j_kp=b->eddyv(i,j,k+1);
	
	visc_ijk=b->visc(i,j,k);
	visc_ip_j_k=b->visc(i+1,j,k);
	visc_i_jm_k=b->visc(i,j-1,k);
	visc_i_jp_k=b->visc(i,j+1,k);
	visc_i_j_km=b->visc(i,j,k-1);
	visc_i_j_kp=b->visc(i,j,k+1);
	
	visc_ddy_p = (vfm*visc_ijk+ev_ijk + vfm*visc_ip_j_k+ev_ip_j_k + vfm*visc_i_jp_k+visc_i_jp_k + vfm*b->visc(i+1,j+1,k)+b->eddyv(i+1,j+1,k))*0.25;
	visc_ddy_m = (vfm*visc_i_jm_k+ev_i_jm_k  +vfm*b->visc(i+1,j-1,k)+b->eddyv(i+1,j-1,k) + vfm*visc_ijk+ev_ijk + vfm*b->visc(i+1,j,k)+b->eddyv(i+1,j,k))*0.25;
	visc_ddz_p = (vfm*visc_ijk+ev_ijk + vfm*visc_ip_j_k+ev_ip_j_k + vfm*visc_i_j_kp+ev_i_j_kp + vfm*b->visc(i+1,j,k+1)+b->eddyv(i+1,j,k+1))*0.25;
	visc_ddz_m = (vfm*b->visc(i,j,k-1)+b->eddyv(i,j,k-1) + vfm*b->visc(i+1,j,k-1)+b->eddyv(i+1,j,k-1) + vfm*visc_ijk+ev_ijk + vfm*visc_ip_j_k+ev_ip_j_k)*0.25;
		
        
	b->M.p[count] =  2.0*(vfm*visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP]*p->DXP[IP])
				   + 2.0*(vfm*visc_ijk+ev_ijk)/(p->DXN[IM1]*p->DXP[IP])
				   + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				   + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				   + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				   + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				   + CPOR1/(alpha*p->dt);
				  
	b->rhsvec.V[count] += ((v(i+1,j,k)-v(i,j,k))*visc_ddy_p - (v(i+1,j-1,k)-v(i,j-1,k))*visc_ddy_m)/(p->DXP[IP]*p->DYN[JP])
						 + ((w(i+1,j,k)-w(i,j,k))*visc_ddz_p - (w(i+1,j,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DXP[IP]*p->DZN[KP])
									
						 + b->M.p[count]*u(i,j,k)*(1.0/p->N54-1.0)
						 + (CPOR1*u(i,j,k))/(alpha*p->dt);
									
	 b->M.p[count] /= p->N54;
	 
	 b->M.s[count] = -2.0*(vfm*visc_ijk+ev_ijk)/(p->DXN[IM1]*p->DXP[IP]);
	 b->M.n[count] = -2.0*(vfm*visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP]*p->DXP[IP]);
	 
	 b->M.e[count] = -visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 b->M.w[count] = -visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
 
	 ++count;
	}
    
    n=0;
	ULOOP
	{
		if(p->flag1[Im1JK]<0)
		{
		b->rhsvec.V[n] -= b->M.s[n]*u(i-1,j,k);
		b->M.s[n] = 0.0;
		}
		
		if(p->flag1[Ip1JK]<0)
		{
		b->rhsvec.V[n] -= b->M.n[n]*u(i+1,j,k);
		b->M.n[n] = 0.0;
		}
		
		if(p->flag1[IJm1K]<0)
		{
		b->rhsvec.V[n] -= b->M.e[n]*u(i,j-1,k);
		b->M.e[n] = 0.0;
		}
		
		if(p->flag1[IJp1K]<0)
		{
		b->rhsvec.V[n] -= b->M.w[n]*u(i,j+1,k);
		b->M.w[n] = 0.0;
		}
		
	++n;
	}
	
	psolv->start(p,a,pgc,u,b->xvec,b->rhsvec,1,gcval_u,p->D29);
    }
	
    pgc->start1(p,u,gcval_u);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
*/
}

void sflow_idiff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{
    SLICELOOP2
    {
        
	b->G(i,j) +=  ((visc+0.5*(b->eddyv(i,j) + b->eddyv(i,j+1)))/(p->dx*p->dx))*
    
                (     (v(i+1,j) - 2.0*v(i,j) + v(i-1,j))
                + 2.0*(v(i,j+1) - 2.0*v(i,j) + v(i,j-1))
                
                + (u(i,j+1)-v(i,j)) - (v(i-1,j+1)-v(i-1,j)));
    }
}