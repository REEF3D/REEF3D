/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"diff_wallghost.h"

void idiff2_FS::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &w_in, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	


    WLOOP
    diff(i,j,k) = w_in(i,j,k);
    
	pgc->start3(p,diff,gcval_w);

    assemble_w(p,a,pgc,w_in,u,v,w,alpha);

	psolv->start(p,a,pgc,diff,a->rhsvec,3);
    
    
	pgc->start3(p,diff,gcval_w);
	
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}

// matrix and rhs of the implicit diffusion of w (into a->M, a->rhsvec, which must be zero)
void idiff2_FS::assemble_w(lexer* p, fdm* a, ghostcell *pgc, field &w_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;
    int n;

	n=0;

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
    
	a->M.p[n] = 2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				  + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				  + CPOR3/(alpha*p->dt);
				  
	a->rhsvec.V[n] +=  ((a->u(i,j,k+1)-u(i,j,k))*visc_ddx_p - (u(i-1,j,k+1)-u(i-1,j,k))*visc_ddx_m)/(p->DZP[KP]*p->DXN[IP])
						+  ((a->v(i,j,k+1)-v(i,j,k))*visc_ddy_p - (v(i,j-1,k+1)-v(i,j-1,k))*visc_ddy_m)/(p->DZP[KP]*p->DYN[JP])
									
						+ (CPOR3*w_in(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[n] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[n] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[n] = -visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[n] = -visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[n] = -2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP]);
	 a->M.t[n] = -2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP]);
     
     if(p->DF3[IJK]<0 && p->D22==1)
     {
        a->M.p[n]  =  1.0;

        a->M.n[n] = 0.0;
        a->M.s[n] = 0.0;

        a->M.w[n] = 0.0;
        a->M.e[n] = 0.0;

        a->M.t[n] = 0.0;
        a->M.b[n] = 0.0;
        
        a->rhsvec.V[n] =  0.0;
     }
    
	 ++n;
	}
    
    if(p->D22==1)
    {
    n=0;
    WLOOP
	{
        if(p->DF3[IJK]>0)
        {
            
		if(p->DF3[Im1JK]<0)  // solid (direct forcing): u = 0 there, Dirichlet
		{
		a->M.s[n] = 0.0;
		}
		else
		if((p->flag3[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0)))
		{
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
		
		if(p->DF3[Ip1JK]<0)  // solid (direct forcing): u = 0 there, Dirichlet
		{
		a->M.n[n] = 0.0;
		}
		else
		if((p->flag3[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0)))
		{
		a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
		
		if(p->DF3[IJm1K]<0)  // solid (direct forcing): u = 0 there, Dirichlet
		{
		a->M.e[n] = 0.0;
		}
		else
		if((p->flag3[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0)))
		{
		a->M.p[n] += a->M.e[n];
		a->M.e[n] = 0.0;
		}
		
		if(p->DF3[IJp1K]<0)  // solid (direct forcing): u = 0 there, Dirichlet
		{
		a->M.w[n] = 0.0;
		}
		else
		if((p->flag3[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0)))
		{
		a->M.p[n] += a->M.w[n];
		a->M.w[n] = 0.0;
		}
		
		if(p->DF3[IJKm1]<0)  // solid (direct forcing): u = 0 there, Dirichlet
		{
		a->M.b[n] = 0.0;
		}
		else
		if((p->flag3[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0)))
		{
		a->M.p[n] += a->M.b[n];
		a->M.b[n] = 0.0;
		}
		
		if(p->DF3[IJKp1]<0)  // solid (direct forcing): u = 0 there, Dirichlet
		{
		a->M.t[n] = 0.0;
		}
		else
		if((p->flag3[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0)))
		{
		a->M.p[n] += a->M.t[n];
		a->M.t[n] = 0.0;
		}
        
        }

	++n;
	}
    }
    
    
    if(p->D22==2)
    {
    // walls: couple the ghost cell to the new values (instead of the old stage field)
    double c1,c2;
    pwall->update(p,pgc,2,gcval_w);

    n=0;
    WLOOP
	{
            
		if(p->flag3[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0))
		{
		if(pwall->coef(p,i,j,k,1,c1,c2) && (c2==0.0 || p->flag3[Ip1JK]>0))
		{
		a->M.p[n] += a->M.s[n]*c1;
		a->M.n[n] += a->M.s[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.s[n]*w(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag3[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0))
		{
		if(pwall->coef(p,i,j,k,4,c1,c2) && (c2==0.0 || p->flag3[Im1JK]>0))
		{
		a->M.p[n] += a->M.n[n]*c1;
		a->M.s[n] += a->M.n[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.n[n]*w(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag3[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0))
		{
		if(pwall->coef(p,i,j,k,3,c1,c2) && (c2==0.0 || p->flag3[IJp1K]>0))
		{
		a->M.p[n] += a->M.e[n]*c1;
		a->M.w[n] += a->M.e[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.e[n]*w(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag3[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0))
		{
		if(pwall->coef(p,i,j,k,2,c1,c2) && (c2==0.0 || p->flag3[IJm1K]>0))
		{
		a->M.p[n] += a->M.w[n]*c1;
		a->M.e[n] += a->M.w[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.w[n]*w(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag3[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0))
		{
		if(pwall->coef(p,i,j,k,5,c1,c2) && (c2==0.0 || p->flag3[IJKp1]>0))
		{
		a->M.p[n] += a->M.b[n]*c1;
		a->M.t[n] += a->M.b[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.b[n]*w(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag3[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0))
		{
		if(pwall->coef(p,i,j,k,6,c1,c2) && (c2==0.0 || p->flag3[IJKm1]>0))
		{
		a->M.p[n] += a->M.t[n]*c1;
		a->M.b[n] += a->M.t[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.t[n]*w(i,j,k+1);
		a->M.t[n] = 0.0;
		}
        
        }

	++n;

    }
}

// explicit evaluation of the same operator at the velocities u, v, w: residual of the
// system assembled with w_in = w, divided by CPOR
void idiff2_FS::apply_w(lexer* p, fdm* a, ghostcell *pgc, field &D, field &u, field &v, field &w)
{
    int q=0;
    WLOOP
    {
    a->rhsvec.V[q]=0.0;
    ++q;
    }

    assemble_w(p,a,pgc,w,u,v,w,1.0);

    q=0;
    WLOOP
    {
    D(i,j,k) = (a->rhsvec.V[q] - (a->M.p[q]*w(i,j,k)
                + a->M.s[q]*w(i-1,j,k)
                + a->M.n[q]*w(i+1,j,k)
                + a->M.e[q]*w(i,j-1,k)
                + a->M.w[q]*w(i,j+1,k)
                + a->M.b[q]*w(i,j,k-1)
                + a->M.t[q]*w(i,j,k+1)))/CPOR3;
    ++q;
    }

    q=0;
    WLOOP
    {
    a->rhsvec.V[q]=0.0;
    ++q;
    }
}

