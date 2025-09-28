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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fluid_update_fsf.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fluid_update_fsf::fluid_update_fsf(lexer *p, fdm* a, ghostcell* pgc) : dx(p->DXM),visc_air(p->W4),visc_water(p->W2),visc_body(p->X44),
                                                                      ro_air(p->W3),ro_water(p->W1),visc_sed(p->W2),ro_sed(p->S22)
{
    gcval_ro=1;
	gcval_visc=1;
}

fluid_update_fsf::~fluid_update_fsf()
{
}

void fluid_update_fsf::start(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
	double H=0.0;
    double H_fb=0.0;
	p->volume1=0.0;
	p->volume2=0.0;
    
    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    
	BASELOOP
	{    
		if(a->phi(i,j,k)>p->psi)
		H=1.0;

		if(a->phi(i,j,k)<-p->psi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=p->psi)
		H=0.5*(1.0 + a->phi(i,j,k)/p->psi + (1.0/PI)*sin((PI*a->phi(i,j,k))/p->psi));


        a->ro(i,j,k)   = ro_water*H +   ro_air*(1.0-H);
        a->visc(i,j,k) = visc_water*H + visc_air*(1.0-H);
            
            if(p->flagsf4[IJK]>0)
            {
            p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
            p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
            }
	}
    

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);
    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
    }
    ++iocheck;
    
    // sediment
    if(p->Q10==1)
    BASELOOP
	{    
		if(a->topo(i,j,k)>p->psi)
		H=1.0;

		if(a->topo(i,j,k)<-p->psi)
		H=0.0;

		if(fabs(a->topo(i,j,k))<=p->psi)
		H=0.5*(1.0 + a->topo(i,j,k)/p->psi + (1.0/PI)*sin((PI*a->topo(i,j,k))/p->psi));


        a->ro(i,j,k)   = a->ro(i,j,k)*H +   ro_sed*(1.0-H);
        a->visc(i,j,k) = a->visc(i,j,k)*H + visc_sed*(1.0-H);
	}
    
    
    pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);
}

int fluid_update_fsf::iocheck;
int fluid_update_fsf::iter;
