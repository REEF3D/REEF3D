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

void idiff2_FS::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &v_in, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
    
    VLOOP
    diff(i,j,k) = v_in(i,j,k);
    
    pgc->start2(p,diff,gcval_v);

    assemble_v(p,a,pgc,v_in,u,v,w,alpha);

	psolv->start(p,a,pgc,diff,a->rhsvec,2);
    
    
	pgc->start2(p,diff,gcval_v);
	
    time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"vdiffiter: "<<p->viter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}

// matrix and rhs of the implicit diffusion of v (into a->M, a->rhsvec, which must be zero)
void idiff2_FS::assemble_v(lexer* p, fdm* a, ghostcell *pgc, field &v_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;
    int n;

	n=0;

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
    
	
	a->M.p[n] = 2.0*(visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP1]*p->DYP[JP])
				  + 2.0*(visc_ijk+ev_ijk)/(p->DYN[JP]*p->DYP[JP])
				  + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
				  + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
				  + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				  + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				  + CPOR2/(alpha*p->dt);
				  
	a->rhsvec.V[n] += ((u(i,j+1,k)-u(i,j,k))*visc_ddx_p - (u(i-1,j+1,k)-u(i-1,j,k))*visc_ddx_m)/(p->DYP[JP]*p->DXN[IP])
						+  ((w(i,j+1,k)-w(i,j,k))*visc_ddz_p - (w(i,j+1,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DYP[JP]*p->DZN[KP])
									
						+ (CPOR2*v_in(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[n] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
	 a->M.n[n] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
	 
	 a->M.e[n] = -2.0*(visc_ijk+ev_ijk)/(p->DYN[JP]*p->DYP[JP]);
	 a->M.w[n] = -2.0*(visc_i_jp_k+ev_i_jp_k)/(p->DYN[JP1]*p->DYP[JP]);
	 
	 a->M.b[n] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[n] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
     
     if(p->DF2[IJK]<0 && p->D22==1)
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
	VLOOP
	{
        if(p->DF2[IJK]>0)
        {
            
		if((p->flag2[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0)) || p->DF2[Im1JK]<0)
		{
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
		
		if((p->flag2[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0)) || p->DF2[Ip1JK]<0)
		{
		a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
		
		if((p->flag2[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0)) || p->DF2[IJm1K]<0)
		{
		a->M.p[n] += a->M.e[n];
		a->M.e[n] = 0.0;
		}
		
		if((p->flag2[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0)) || p->DF2[IJp1K]<0)
		{
		a->M.p[n] += a->M.w[n];
		a->M.w[n] = 0.0;
		}
		
		if((p->flag2[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0)) || p->DF2[IJKm1]<0)
		{
		a->M.p[n] += a->M.b[n];
		a->M.b[n] = 0.0;
		}
		
		if((p->flag2[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0)) || p->DF2[IJKp1]<0)
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
    pwall->update(p,pgc,1,gcval_v);

    n=0;
	VLOOP
	{
            
		if(p->flag2[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0))
		{
		if(pwall->coef(p,i,j,k,1,c1,c2) && (c2==0.0 || p->flag2[Ip1JK]>0))
		{
		a->M.p[n] += a->M.s[n]*c1;
		a->M.n[n] += a->M.s[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.s[n]*v(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag2[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0))
		{
		if(pwall->coef(p,i,j,k,4,c1,c2) && (c2==0.0 || p->flag2[Im1JK]>0))
		{
		a->M.p[n] += a->M.n[n]*c1;
		a->M.s[n] += a->M.n[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.n[n]*v(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag2[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0))
		{
		if(pwall->coef(p,i,j,k,3,c1,c2) && (c2==0.0 || p->flag2[IJp1K]>0))
		{
		a->M.p[n] += a->M.e[n]*c1;
		a->M.w[n] += a->M.e[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.e[n]*v(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag2[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0))
		{
		if(pwall->coef(p,i,j,k,2,c1,c2) && (c2==0.0 || p->flag2[IJm1K]>0))
		{
		a->M.p[n] += a->M.w[n]*c1;
		a->M.e[n] += a->M.w[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.w[n]*v(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag2[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0))
		{
		if(pwall->coef(p,i,j,k,5,c1,c2) && (c2==0.0 || p->flag2[IJKp1]>0))
		{
		a->M.p[n] += a->M.b[n]*c1;
		a->M.t[n] += a->M.b[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.b[n]*v(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag2[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0))
		{
		if(pwall->coef(p,i,j,k,6,c1,c2) && (c2==0.0 || p->flag2[IJKm1]>0))
		{
		a->M.p[n] += a->M.t[n]*c1;
		a->M.b[n] += a->M.t[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.t[n]*v(i,j,k+1);
		a->M.t[n] = 0.0;
		}
        

	++n;
	}
    }
}

// explicit evaluation of the same operator at the velocities u, v, w: residual of the
// system assembled with v_in = v, divided by CPOR
void idiff2_FS::apply_v(lexer* p, fdm* a, ghostcell *pgc, field &D, field &u, field &v, field &w)
{
    int q=0;
    VLOOP
    {
    a->rhsvec.V[q]=0.0;
    ++q;
    }

    assemble_v(p,a,pgc,v,u,v,w,1.0);

    q=0;
    VLOOP
    {
    D(i,j,k) = (a->rhsvec.V[q] - (a->M.p[q]*v(i,j,k)
                + a->M.s[q]*v(i-1,j,k)
                + a->M.n[q]*v(i+1,j,k)
                + a->M.e[q]*v(i,j-1,k)
                + a->M.w[q]*v(i,j+1,k)
                + a->M.b[q]*v(i,j,k-1)
                + a->M.t[q]*v(i,j,k+1)))/CPOR2;
    ++q;
    }

    q=0;
    VLOOP
    {
    a->rhsvec.V[q]=0.0;
    ++q;
    }
}


