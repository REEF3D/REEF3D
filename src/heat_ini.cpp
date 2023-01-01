/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"heat_print.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"

void heat_print::heat_ini(lexer* p, fdm *a, ghostcell* pgc,heat *pheat)
{
	if(p->H10>0 && p->W90==0 && p->H3==1)
	pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);
    
    if(p->H10>0 && p->W90==0 && p->H3==2)
	pupdate = new fluid_update_fsf_heat_Bouss(p,a,pgc,pheat);

double dx=p->DXM;
double r;


    LOOP
	T(i,j,k)=p->H50_2;

    double psi=1.0e-20;

	LOOP
	if(p->XP[IP]>p->H51-psi && p->XP[IP]<p->H54+psi
	&& p->YP[JP]>p->H52-psi && p->YP[JP]<p->H55+psi
	&& p->ZP[KP]>p->H53-psi && p->ZP[KP]<p->H56+psi)
	T(i,j,k)=p->H50_1;



    if(p->H57_1>0||p->H57_2>0||p->H57_3>0||p->H57_4>0)
    {
        LOOP
        if(p->H57_1*p->XP[IP]+ p->H57_2*p->YP[JP] + p->H57_3*p->ZP[KP] < p->H57_4)
        T(i,j,k)=p->H50_1;
    }
    

    if(p->H58_4>0.0)
    {

        LOOP
        {
        r = sqrt( pow(p->XP[IP]-p->H58_1,2.0)+pow(p->YP[JP]-p->H58_2,2.0)+pow(p->ZP[KP]-p->H58_3,2.0));
        if(r<=p->H58_4)
        T(i,j,k)=p->H50_1;
        }
    }

    pgc->start4(p,T,80); 
    pgc->start4(p,T,80);
    
    pupdate->start(p,a,pgc);
    pgc->start4(p,a->ro,1);

}
