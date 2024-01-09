/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"solver2D.h"

void sediment_exner::non_equillibrium_solve(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    double rhosed=p->S22;
    double rhowat=p->W1;
    double g=9.81;
    double d50=p->S20;
    double visc=p->W2;
    double kappa=0.4;
    double ks=p->S21*d50;
    double Rstar=(rhosed-rhowat)/rhowat;
    double Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
    double Ti;
    double time,starttime;
    
    
    SLICELOOP4
    {
    Ti=MAX((s->shearvel_eff(i,j)*s->shearvel_eff(i,j)-s->shearvel_crit(i,j)*s->shearvel_crit(i,j))/(s->shearvel_crit(i,j)*s->shearvel_crit(i,j)),0.0);
        
    //Ls = 3.0*d50*pow(Ds,0.6)*pow(Ti,0.9);
    
    Ls = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j), 0.0)*d50;
    
    //Ls = p->dtsed/p->DXM*sqrt(pow(0.5*(s->P(i,j)+s->P(i+1,j)),2.0) +  pow(0.5*(s->Q(i,j)+s->Q(i,j+1)),2.0));
    
    Ls = MAX(Ls,0.0);
    Ls = MIN(Ls,1.0);

    s->qb(i,j) =  s->qbe(i,j) - Ls*(dqx0(i,j) + dqy0(i,j));
    }
    pgc->gcsl_start4(p,s->qb,1);
    
    
// Implicit Solution

/*
starttime=pgc->timer();

    n=0;
    SLICELOOP4
	{
    Ls = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j), 1.0e-6)*d50;
    
	M.p[n]  =  1.0;     
               

   	M.n[n] = Ls/(p->DXP[IP]+p->DXP[IM1]);
	M.s[n] = -Ls/(p->DXP[IP]+p->DXP[IM1]);

	M.w[n] = Ls/(p->DYP[JP]+p->DYP[JM1])*p->y_dir;
	M.e[n] = -Ls/(p->DYP[JP]+p->DYP[JM1])*p->y_dir;
    
    rhsvec.V[n] = s->qbe(i,j);
    s->qb(i,j)  = s->qbe(i,j);
    
	++n;
	}
	
	
    n=0;
	SLICELOOP4
	{
		if(p->flagslice4[Im1J]<0)
		{
		rhsvec.V[n] -= M.s[n]*s->qb(i-1,j);
		M.s[n] = 0.0;
		}
		
		if(p->flagslice4[Ip1J]<0)
		{
		rhsvec.V[n] -= M.n[n]*s->qb(i+1,j);
		M.n[n] = 0.0;
		}
		
		if(p->flagslice4[IJm1]<0)
		{
		rhsvec.V[n] -= M.e[n]*s->qb(i,j-1);
		M.e[n] = 0.0;
		}
		
		if(p->flagslice4[IJp1]<0)
		{
		rhsvec.V[n] -= M.w[n]*s->qb(i,j+1);
		M.w[n] = 0.0;
		}
		
	++n;
	}
    
    psolv->start(p,pgc,qb,M,xvec,rhsvec,4);
    
    pgc->gcsl_start4(p,qb,1);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"qb_iter: "<<p->uiter<<"  qb_time: "<<setprecision(3)<<time<<endl;*/
}