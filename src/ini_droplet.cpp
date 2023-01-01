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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::droplet_ini(lexer* p, fdm *a, ghostcell* pgc)
{

double dx=p->DXM;
double r;
double diff = fabs(p->I58_2-p->F58_4)>1.0e-9?p->I58_2-p->F58_4 : 1.0e20;

if(p->F58_4>0.0)
{

	WLOOP
	{
    r = sqrt( pow(p->XP[IP]-p->F58_1,2.0)+pow(p->YP[JP]-p->F58_2,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));

         if(r<=p->F58_4)
        {
        a->w(i,j,k)=p->I58_1;
        }

        if(r<=p->I58_2 && r>p->F58_4)
        {
        a->w(i,j,k)=p->I58_1 - p->I58_1*(r-p->F58_4)/(diff);
        }
	}
}
	pgc->start3(p,a->w,12);



    LOOP
	{
        epsi = (1.6/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
        
		if(a->phi(i,j,k)>=0)
		H=1.0;

		if(a->phi(i,j,k)<0)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}

}


