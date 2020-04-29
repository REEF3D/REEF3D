/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"fluid_update_fsf.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fluid_update_fsf::fluid_update_fsf(lexer *p, fdm* a, ghostcell* pgc) : dx(p->dx),visc_air(p->W4),visc_water(p->W2),
                                                                      ro_air(p->W3),ro_water(p->W1)
{
    gcval_ro=1;
	gcval_visc=1;
}

fluid_update_fsf::~fluid_update_fsf()
{
}

void fluid_update_fsf::start(lexer *p, fdm* a, ghostcell* pgc)
{
	double H=0.0;
	p->volume1=0.0;
	p->volume2=0.0;
    
    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    
	LOOP
	{
  
       epsi = p->F45*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		a->ro(i,j,k)=     ro_water*H +   ro_air*(1.0-H);
		a->visc(i,j,k) = visc_water*H + visc_air*(1.0-H);
        
        a->visctot(i,j,k) = a->visc(i,j,k) + a->eddyv(i,j,k);

		p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
		p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}
    
	pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);
    pgc->start4(p,a->visctot,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);
    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
    }
    ++iocheck;

}

int fluid_update_fsf::iocheck;
int fluid_update_fsf::iter;
