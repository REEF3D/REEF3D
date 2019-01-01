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

#include"fluid_update_fsf_entrain.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"concentration.h"

fluid_update_fsf_entrain::fluid_update_fsf_entrain(lexer *p, fdm* a, ghostcell* pgc, concentration *&ppconcentration) : dx(p->dx)
{
    gcval_ro=1;
	gcval_visc=1;

	visc_air = p->W4;
	visc_water = p->W2;
	ro_air = p->W3;
	ro_water = p->W1;
	
	pconcentration = ppconcentration;
}

fluid_update_fsf_entrain::~fluid_update_fsf_entrain()
{
}

void fluid_update_fsf_entrain::start(lexer *p, fdm* a, ghostcell* pgc)
{
	double H=0.0;
	double conc;
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
		
		conc=pconcentration->val(i,j,k);

		a->ro(i,j,k)=      (ro_water-conc*ro_water + conc*ro_air)*H +   (ro_air)*(1.0-H);
		
		a->visc(i,j,k)=    (visc_water)*H + (visc_air)*(1.0-H);

		p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
		p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}

	pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);

    if(p->mpirank==0 && iocheck==0 && (innercounter==p->N50-1 || p->N52==0))
    {
	cout<<"Volume_entr 1: "<<p->volume1<<endl;
	cout<<"Volume_entr 2: "<<p->volume2<<endl;
    }
    ++iocheck;

}

void fluid_update_fsf_entrain::start3(lexer *p, fdm* a, ghostcell* pgc, field &ls1, field &ls2)
{
}

int fluid_update_fsf_entrain::iocheck;
int fluid_update_fsf_entrain::iter;
