/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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
--------------------------------------------------------------------*/

#include<sys/stat.h>
#include"mooring_Catenary.h"
#include"lexer.h"
#include"ghostcell.h"


void mooring_Catenary::initialize(lexer *p, ghostcell *pgc)
{
	double rho_f = 1000.0;
	
	rho_c = p->X311_rho_c[line];
	w = p->X311_w[line]*9.81*(rho_c - rho_f)/rho_c;
	L = p->X311_l[line];
	H = p->X311_H[line];
	EA = p->X311_EA[line];
	
	xs = p->X311_xs[line];
	ys = p->X311_ys[line];
	zs = p->X311_zs[line];
	
	p->Darray(x,H); 
	p->Darray(y,H);
	p->Darray(z,H); 
	p->Darray(T,H);
	p->Darray(B,2);
	p->Darray(F,2);
	p->Darray(A,2,2);

	if(p->mpirank==0)
	{
		char str[1000];
		sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;	
	}
	printtime = 0.0;
    
    // Initialise breaking
    broken = false;
    curr_time = 0.0;
    breakTension = p->X314 > 0 ? p->X314_T[line]: 0.0;
    breakTime = p->X315 > 0 ? p->X315_t[line]: 0.0;
}
