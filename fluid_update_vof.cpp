/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"fluid_update_vof.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fluid_update_vof::fluid_update_vof(lexer *p, fdm* a, ghostcell* pgc) : dx(p->dx),
												visc_air(p->W4),visc_water(p->W2),ro_air(p->W3),ro_water(p->W1)
{
    gcval_ro=1;
	gcval_visc=1;
}

fluid_update_vof::~fluid_update_vof()
{
}

void fluid_update_vof::start(lexer *p, fdm* a, ghostcell* pgc)
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
        
		H=a->vof(i,j,k);

		H=MAX(H,0.0);
		H=MIN(H,1.0);

		a->ro(i,j,k)=      ro_water*H +   ro_air*(1.0-H);
		a->visc(i,j,k)= visc_water*H + visc_air*(1.0-H);

		p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
		p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}

	pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);
    
    if(p->mpirank==0 && iocheck==0 && (innercounter==p->N50-1 || p->N52==0) && (p->count%p->P12==0))
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
    }
    ++iocheck;

}

void fluid_update_vof::start3(lexer *p, fdm* a, ghostcell* pgc, field &ls1, field &ls2)
{
}

int fluid_update_vof::iocheck;
int fluid_update_vof::iter;
