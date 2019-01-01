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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::morison(lexer* p, fdm *a, ghostcell *pgc)
{
    //cd = morison_cd(p,a,pgc,Re,Kc);
    //cm = morison_cm(p,a,pgc,Re,Kc);
		
	// cout<<"Cm: "<<cm<<endl;
	
	cm=p->P89_cm;
	cd=p->P89_cd;
	
	
	F_morison = p->W1*cm*Vcyl*uinf_dt
               + 0.5*p->W1*cd*Ai*uinf*fabs(uinf);

    Fmi = p->W1*cm*Vcyl*vinf_dt;
    Fmd = 0.5*p->W1*cd*Ai*uinf*fabs(uinf);

    F_morison_norm = F_morison/(p->W1*fabs(p->W22)*p->wH*p->wH*p->P82_y);
	
	
	F_morison_rect = p->W1*cm*Vrect*uinf_dt
               + 0.5*p->W1*cd*Ai*uinf*fabs(uinf);
}


