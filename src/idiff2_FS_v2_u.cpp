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

#include"idiff2_FS_v2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

idiff2_FS_v2::idiff2_FS_v2(lexer* p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    if(p->B21==1)
    {
    gcval_udiff=10;
	gcval_vdiff=11;
	gcval_wdiff=12;
    }
    
    if(p->B21==2)
    {
    gcval_udiff=117;
	gcval_vdiff=118;
	gcval_wdiff=119;
    }
    
    if(p->B21==3)
    {
    gcval_udiff=110;
	gcval_vdiff=111;
	gcval_wdiff=112;
    }
}

idiff2_FS_v2::~idiff2_FS_v2()
{
}

void idiff2_FS_v2::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	double visc;
    
    pgc->start1(p,u,gcval_udiff);
	pgc->start2(p,v,gcval_vdiff);
	pgc->start3(p,w,gcval_wdiff);

    count=0;

    
	count=0;
    if(p->i_dir==1)
    {
	ULOOP
	{
	visc = 0.5*(a->visc(i,j,k) + a->visc(i+1,j,k)) + 0.5*(a->eddyv(i,j,k) + a->eddyv(i+1,j,k));
    
	a->M.p[count] =  2.0*visc/(p->DXN[IP1]*p->DXP[IP])
				   + 2.0*visc/(p->DXN[IP]*p->DXP[IP])
				   + visc/(p->DYP[JP]*p->DYN[JP])
				   + visc/(p->DYP[JM1]*p->DYN[JP])
				   + visc/(p->DZP[KP]*p->DZN[KP])
				   + visc/(p->DZP[KM1]*p->DZN[KP])
				   + CPOR1/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((v(i+1,j,k)-v(i,j,k))*visc - (v(i+1,j-1,k)-v(i,j-1,k))*visc)/(p->DXP[IP]*p->DYN[JP])
						 + ((w(i+1,j,k)-w(i,j,k))*visc - (w(i+1,j,k-1)-w(i,j,k-1))*visc)/(p->DXP[IP]*p->DZN[KP])

						 + (CPOR1*u(i,j,k))/(alpha*p->dt);
                         
	 
	 a->M.s[count] = -2.0*visc/(p->DXN[IP]*p->DXP[IP]);
	 a->M.n[count] = -2.0*visc/(p->DXN[IP1]*p->DXP[IP]);
	 
	 a->M.e[count] = -visc/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] = -visc/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] = -visc/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc/(p->DZP[KP]*p->DZN[KP]);
	 
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
		
		if(p->flag1[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*u(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag1[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*u(i,j+1,k);
		a->M.w[n] = 0.0;
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
	
	psolv->start(p,a,pgc,u,a->rhsvec,1);
    }
	
    pgc->start1(p,u,gcval_u);
    
    pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
}

