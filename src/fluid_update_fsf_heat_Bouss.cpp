/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"fluid_update_fsf_heat_Bouss.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heat.h"

fluid_update_fsf_heat_Bouss::fluid_update_fsf_heat_Bouss(lexer *p, fdm* a, ghostcell* pgc, heat *&ppheat) : dx(p->DXM)
{
    gcval_ro=1;
	gcval_visc=1;

	visc_2 = p->W4;
	visc_1 = p->W2;
	ro_2 = p->W3;
	ro_1 = p->W1;
	alpha_air = p->H2;
	alpha_water = p->H1;
    
    if(p->H9==1)
    {
    T0_1 = p->H50_1 + 273.0;
    T0_2 = p->H50_2 + 273.0;
    }
    
    if(p->H9==2)
    {
    T0_1 = p->H50_2 + 273.0;
    T0_2 = p->H50_1 + 273.0;
    }
	
	pheat = ppheat;
    
}

fluid_update_fsf_heat_Bouss::~fluid_update_fsf_heat_Bouss()
{
}

void fluid_update_fsf_heat_Bouss::start(lexer *p, fdm* a, ghostcell* pgc)
{
    
	double H=0.0;
	double temp;
	p->volume1=0.0;
	p->volume2=0.0;

    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    
    if(p->j_dir==0)        
    epsi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);

   //
	LOOP
	{
        temp = pheat->val(i,j,k) + 273.0;
        
        if(p->H4==0)
        {
            if(p->H9==1)
            {
            ro_1 = p->W1 - p->W1*(temp - T0_1)/T0_1;
            ro_2 = p->W3 - p->W3*(temp - T0_2)/T0_2;

            visc_1 = p->W2;
            visc_2 = p->W4;
            }
            
            if(p->H9==2)
            {
            ro_1 = p->W3 - p->W3*(temp - T0_2)/T0_2;
            ro_2 = p->W1 - p->W1*(temp - T0_1)/T0_1;

            visc_1 = p->W4;
            visc_2 = p->W2;
            }
        }
        
        if(p->H4==1)
        {
            if(p->H9==1)
            {
            ro_1 = p->W1 - p->W1*p->H4_beta1*(temp - T0_1);
            ro_2 = p->W3 - p->W3*p->H4_beta2*(temp - T0_2);

            visc_1 = p->W2;
            visc_2 = p->W4;
            }
            
            if(p->H9==2)
            {
            ro_1 = p->W3 - p->W3*p->H4_beta2*(temp - T0_2);
            ro_2 = p->W1 - p->W1*p->H4_beta1*(temp - T0_1);

            visc_1 = p->W4;
            visc_2 = p->W2;
            }
        }

		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		a->ro(i,j,k)=     ro_1*H +   ro_2*(1.0-H);
		a->visc(i,j,k)= visc_1*H + visc_2*(1.0-H);

		p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
		p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}

	pgc->start4(p,a->ro,gcval_ro);
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

int fluid_update_fsf_heat_Bouss::iocheck;
int fluid_update_fsf_heat_Bouss::iter;
