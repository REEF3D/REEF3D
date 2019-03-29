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
#include"ghostcell.h"
#include"fdm2D.h"
#include"solver2D.h"

sflow_idiff::sflow_idiff(lexer* p)
{
    gcval_u = 10;
    gcval_v = 11;
}

sflow_idiff::~sflow_idiff()
{
}

void sflow_idiff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{
    starttime=pgc->timer();
   
    SLICELOOP1
    {
	visc = p->W2 + 0.5*(b->eddyv(i,j) + b->eddyv(i+1,j));

        
	b->M.p[count] =  6.0*visc/(p->dx*p->dx);
                   
				   + 1.0/(alpha*p->dt);
				  
	b->rhsvec.V[count] += (visc/(p->dx*p->dx))*((v(i+1,j)-v(i,j)) - (v(i+1,j-1)-v(i,j-1)))
									
						 + b->M.p[count]*u(i,j)*(1.0/p->N54-1.0)
						 + (u(i,j))/(alpha*p->dt);
									
	 b->M.p[count] /= p->N54;
	 
	 b->M.s[count] = -2.0*visc/(p->dx*p->dx);
	 b->M.n[count] = -2.0*visc/(p->dx*p->dx);
	 
	 b->M.e[count] = -visc/(p->dx*p->dx);
	 b->M.w[count] = -visc/(p->dx*p->dx);
 
	 ++count;
	}
    
    
	psolv->start(p,b,pgc,u,b->xvec,b->rhsvec,1,gcval_u,p->D29);
    
    pgc->gcsl_start1(p,u,gcval_u);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
}

void sflow_idiff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{
    starttime=pgc->timer();
   
    SLICELOOP2
    {
	visc = p->W2 + 0.5*(b->eddyv(i,j) + b->eddyv(i,j+1));

        
	b->M.p[count] =  6.0*visc/(p->dx*p->dx);
                   
				   + 1.0/(alpha*p->dt);
				  
	b->rhsvec.V[count] += (visc/(p->dx*p->dx))*((u(i,j+1)-u(i,j)) - (v(i-1,j+1)-v(i-1,j)))
									
						 + b->M.p[count]*v(i,j)*(1.0/p->N54-1.0)
						 + (v(i,j))/(alpha*p->dt);
									
	 b->M.p[count] /= p->N54;
	 
	 b->M.s[count] = -visc/(p->dx*p->dx);
	 b->M.n[count] = -visc/(p->dx*p->dx);
	 
	 b->M.e[count] = -2.0*visc/(p->dx*p->dx);
	 b->M.w[count] = -2.0*visc/(p->dx*p->dx);
 
	 ++count;
	}
    
    
	psolv->start(p,b,pgc,v,b->xvec,b->rhsvec,2,gcval_v,p->D29);
    
    pgc->gcsl_start2(p,v,gcval_v);
    
	time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
	cout<<"vdiffiter: "<<p->uiter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}