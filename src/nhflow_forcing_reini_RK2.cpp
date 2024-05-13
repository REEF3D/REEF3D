/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_reinidisc_fsf.h"

void nhflow_forcing::reini_RK2(lexer* p, fdm_nhf* d, ghostcell* pgc, double *F)
{	
    LOOP
	dt[IJK] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]*p->sigz[IJ]);

	reiniter=5;
	
	
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<endl<<"initializing reini forcing..."<<endl<<endl;
	reiniter=5;
	}

    for(int q=0;q<reiniter;++q)
	{
        // Step 1
		prdisc->start(p,pgc,F,L);

		LOOP
		FRK1[IJK] = F[IJK] + dt[IJK]*L[IJK];

         pgc->start5V(p,F,1);
        
        
        // Step 2
		prdisc->start(p,pgc,F,L);

		LOOP
		F[IJK] = 0.5*F[IJK] + 0.5*FRK1[IJK] + 0.5*dt[IJK]*L[IJK];

        pgc->start5V(p,F,1);
	}
}





