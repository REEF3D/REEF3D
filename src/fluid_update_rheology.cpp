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

#include"fluid_update_rheology.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"rheology_f.h"

fluid_update_rheology::fluid_update_rheology(lexer *p, fdm* a) : dx(p->DXM),
												visc2(p->W4),ro2(p->W3),ro1(p->W1)
{
    gcval_ro=1;
	gcval_visc=1;
	
	prheo = new rheology_f(p,a);	
}

fluid_update_rheology::~fluid_update_rheology()
{
}

void fluid_update_rheology::start(lexer *p, fdm* a, ghostcell* pgc)
{
	double H=0.0;
	double Hro=0.0;
	p->volume1=0.0;
	p->volume2=0.0;
    
    if(p->count>iter)
    iocheck=0;
	iter=p->count;
	
	visc1 = p->W2;
    
    if(p->j_dir==0)        
    epsi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);


	// density
	LOOP
	{
		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		a->ro(i,j,k)=      ro1*H +   ro2*(1.0-H);

        p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
        p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}
	pgc->start4(p,a->ro,gcval_ro);
    
	// viscosity
	LOOP
	{  
		if(a->phi(i,j,k)>epsi)
		{
		H=1.0;
		visc1 = prheo->viscosity(p,a,pgc);
		}

		if(a->phi(i,j,k)<-epsi)
        {
		H=0.0;
        visc1 = 0.0;
        }

		if(fabs(a->phi(i,j,k))<=epsi)
		{
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));
		visc1 = prheo->viscosity(p,a,pgc);
		}

		a->visc(i,j,k) =    visc1*H + visc2*(1.0-H);
	}
	
	pgc->start4(p,a->visc,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);
    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
    }
    ++iocheck;
}

int fluid_update_rheology::iocheck;
int fluid_update_rheology::iter;

