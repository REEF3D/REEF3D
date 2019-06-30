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

#include"fnpf_sg_fsfbc.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"solver2D.h"

void fnpf_sg_fsfbc::damping(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f, int gcval, double alpha)
{
    double starttime=pgc->timer();
    
    n=0;
    SLICELOOP4
    {
	visc = c->vb(i,j);

        
	c->M.p[n] =   visc/(p->DXP[IM1]*p->DXN[IP])
                +visc/(p->DXP[IP]*p->DXN[IP])
                +visc/(p->DYP[JM1]*p->DYN[JP])
                +visc/(p->DYP[JP]*p->DYN[JP])
                   
				   + 1.0/(alpha*p->dt);
    
	c->rhsvec.V[n] =    c->M.p[n]*f(i,j)*(1.0/p->N54-1.0)
                         
						 + (f(i,j))/(alpha*p->dt);

	 c->M.p[n] /= p->N54;
	 
	 c->M.s[n] = -visc/(p->DXP[IM1]*p->DXN[IP]);
	 c->M.n[n] = -visc/(p->DXP[IP]*p->DXN[IP]);
	 
	 c->M.e[n] = -visc/(p->DYP[JM1]*p->DYN[JP]);
	 c->M.w[n] = -visc/(p->DYP[JP]*p->DYN[JP]);
 
	 ++n;
	}
    
    
    n=0;
    SLICELOOP4
    {
        if(p->flagslice1[Im1J]<0)
		{
		c->rhsvec.V[n] -= c->M.s[n]*f(i-1,j);
		c->M.s[n] = 0.0;
		}
		
		if(p->flagslice1[Ip1J]<0)
		{
		c->rhsvec.V[n] -= c->M.n[n]*f(i+1,j);
		c->M.n[n] = 0.0;
		}
		
		if(p->flagslice1[IJm1]<0)
		{
		c->rhsvec.V[n] -= c->M.e[n]*f(i,j-1);
		c->M.e[n] = 0.0;
		}
		
		if(p->flagslice1[IJp1]<0)
		{
		c->rhsvec.V[n] -= c->M.w[n]*f(i,j+1);
		c->M.w[n] = 0.0;
		}
 
	++n;
	}
    
	//psolv->start(p,c,pgc,f,c->xvec,c->rhsvec,1,50,p->D29);
    
    pgc->gcsl_start4(p,f,gcval);
    
	//time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
	cout<<"fsfbc_damping: "<<p->uiter<<"  fsfbc_damping_time: "<<setprecision(3)<<time<<endl;
    
    
}