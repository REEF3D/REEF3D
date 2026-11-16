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

#include"idiff2_FS_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"
#include"diff_wallghost.h"
#include"diff_walllaw.h"
#include"diff_dynamic.h"

idiff2_FS_2D::idiff2_FS_2D(lexer* p)
{
    pwall = new diff_wallghost(p);
    pwl = new diff_walllaw;
    pdyn = new diff_dynamic(p);

	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

idiff2_FS_2D::~idiff2_FS_2D()
{
    delete pwall;
    delete pwl;
    delete pdyn;
}

void idiff2_FS_2D::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &u_in, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
    
    ULOOP
    diff(i,j,k) = u_in(i,j,k);
    
    pgc->start1(p,diff,gcval_u);

    assemble_u(p,a,pgc,u_in,u,v,w,alpha);

	psolv->start(p,a,pgc,diff,a->rhsvec,1);

	
    pgc->start1(p,diff,gcval_u);

    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
}

// matrix and rhs of the implicit diffusion of u (into a->M, a->rhsvec, which must be zero)
void idiff2_FS_2D::assemble_u(lexer* p, fdm* a, ghostcell *pgc, field &u_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddz_p,visc_ddz_m;


	n=0;

	// D 24 > 0: the dynamic form, mu = rho (nu + nu_t) and 1/rho at the velocity node (diff_dynamic.h)
    pdyn->update(p,a);
    field &vi = pdyn->visc(a);
    field &ed = pdyn->eddyv(a);
    double ir;

    ULOOP
	{
	ev_ijk=ed(i,j,k);
	ev_ip_j_k=ed(i+1,j,k);
	ev_i_j_km=ed(i,j,k-1);
	ev_i_j_kp=ed(i,j,k+1);
	
	visc_ijk=vi(i,j,k);
	visc_ip_j_k=vi(i+1,j,k);
	visc_i_j_km=vi(i,j,k-1);
	visc_i_j_kp=vi(i,j,k+1);
	
	visc_ddz_p = (visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k + visc_i_j_kp+ev_i_j_kp + vi(i+1,j,k+1)+ed(i+1,j,k+1))*0.25;
	visc_ddz_m = (vi(i,j,k-1)+ed(i,j,k-1) + vi(i+1,j,k-1)+ed(i+1,j,k-1) + visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k)*0.25;

    if(pdyn->mode==2)
    {
    visc_ddz_p = pdyn->corner(visc_ijk+ev_ijk, visc_ip_j_k+ev_ip_j_k, visc_i_j_kp+ev_i_j_kp, vi(i+1,j,k+1)+ed(i+1,j,k+1));
    visc_ddz_m = pdyn->corner(vi(i,j,k-1)+ed(i,j,k-1), vi(i+1,j,k-1)+ed(i+1,j,k-1), visc_ijk+ev_ijk, visc_ip_j_k+ev_ip_j_k);
    }

    ir = pdyn->irho(a,0,i,j,k);

		
        
	a->M.p[n] = (2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP])
				   + 2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP])
				   + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				   + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]))*ir
				   + CPOR1/(alpha*p->dt);
				  
	a->rhsvec.V[n] += (((w(i+1,j,k)-w(i,j,k))*visc_ddz_p - (w(i+1,j,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DXP[IP]*p->DZN[KP]))*ir

						 + (CPOR1*u_in(i,j,k))/(alpha*p->dt);

	 
	 a->M.s[n] = -2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP])*ir;
	 a->M.n[n] = -2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP])*ir;
	 
	 a->M.b[n] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])*ir;
	 a->M.t[n] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP])*ir;
     
     if(p->DF1[IJK]<0 && p->D22==1)
     {
        a->M.p[n]  =  1.0;

        a->M.n[n] = 0.0;
        a->M.s[n] = 0.0;

        a->M.t[n] = 0.0;
        a->M.b[n] = 0.0;
        
        a->rhsvec.V[n] =  0.0;
     }
	 
	 ++n;
	}
    
    // wall-law faces (B 10): zero gradient, the wall law of bcmom is the only wall stress there (diff_walllaw.h);
    // D 22 1 cuts these faces below
    if(p->D22!=1)
    {
    pwl->update(p,0);
    pwl->rows(p,a->M,0);
    }
    
    if(p->D22==1)
    {
    n=0;
	ULOOP
	{
        if(p->DF1[IJK]>0)
        {
            
		if((p->flag1[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0)) || p->DF1[Im1JK]<0)
		{
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
		
		if((p->flag1[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0)) || p->DF1[Ip1JK]<0)
		{
		a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
		
		if((p->flag1[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0)) || p->DF1[IJKm1]<0)
		{
		a->M.p[n] += a->M.b[n];
		a->M.b[n] = 0.0;
		}
		
		if((p->flag1[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0)) || p->DF1[IJKp1]<0)
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
    pwall->update(p,pgc,0,gcval_u);

    n=0;
	ULOOP
	{
		if(p->flag1[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0))
		{
		if(pwall->coef(p,i,j,k,1,c1,c2) && (c2==0.0 || p->flag1[Ip1JK]>0))
		{
		a->M.p[n] += a->M.s[n]*c1;
		a->M.n[n] += a->M.s[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.s[n]*u(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag1[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0))
		{
		if(pwall->coef(p,i,j,k,4,c1,c2) && (c2==0.0 || p->flag1[Im1JK]>0))
		{
		a->M.p[n] += a->M.n[n]*c1;
		a->M.s[n] += a->M.n[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.n[n]*u(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag1[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0))
		{
		if(pwall->coef(p,i,j,k,5,c1,c2) && (c2==0.0 || p->flag1[IJKp1]>0))
		{
		a->M.p[n] += a->M.b[n]*c1;
		a->M.t[n] += a->M.b[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.b[n]*u(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag1[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0))
		{
		if(pwall->coef(p,i,j,k,6,c1,c2) && (c2==0.0 || p->flag1[IJKm1]>0))
		{
		a->M.p[n] += a->M.t[n]*c1;
		a->M.b[n] += a->M.t[n]*c2;
		}
		else
		a->rhsvec.V[n] -= a->M.t[n]*u(i,j,k+1);
		a->M.t[n] = 0.0;
		}
        

	++n;
	}
    }
}

// explicit evaluation of the same operator at the velocities u, v, w: residual of the
// system assembled with u_in = u, divided by CPOR
void idiff2_FS_2D::apply_u(lexer* p, fdm* a, ghostcell *pgc, field &D, field &u, field &v, field &w)
{
    int q=0;
    ULOOP
    {
    a->rhsvec.V[q]=0.0;
    ++q;
    }

    assemble_u(p,a,pgc,u,u,v,w,1.0);

    q=0;
    ULOOP
    {
    D(i,j,k) = (a->rhsvec.V[q] - (a->M.p[q]*u(i,j,k)
                + a->M.s[q]*u(i-1,j,k)
                + a->M.n[q]*u(i+1,j,k)
                + a->M.b[q]*u(i,j,k-1)
                + a->M.t[q]*u(i,j,k+1)))/CPOR1;
    ++q;
    }

    q=0;
    ULOOP
    {
    a->rhsvec.V[q]=0.0;
    ++q;
    }
}
