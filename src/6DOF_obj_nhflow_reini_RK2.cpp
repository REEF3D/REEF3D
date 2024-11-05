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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_reinidisc_fsf.h"

void sixdof_obj::nhflow_reini_RK2(lexer* p, fdm_nhf* d, ghostcell* pgc, double *F)
{	
    if(p->j_dir==0)
    LOOP
	DTT[IJK] = 0.5*MIN(p->DXP[IP],p->DZP[KP]/p->sigz[IJ]);
    
    if(p->j_dir==1)
    LOOP
	DTT[IJK] = 0.5*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]/p->sigz[IJ]);

	reiniter=5;
	
	if(p->count==0 && p->mpirank==0)
	cout<<endl<<"initializing reini forcing..."<<endl<<endl;


    for(int q=0;q<reiniter;++q)
	{
        // Step 1
		pnhfrdisc->start(p,pgc,F,LL);

		LOOP
		FRK1[IJK] = F[IJK] + DTT[IJK]*LL[IJK];

         pgc->start5V(p,FRK1,1);
        
        
        // Step 2
		pnhfrdisc->start(p,pgc,FRK1,LL);

		LOOP
		F[IJK] = 0.5*F[IJK] + 0.5*FRK1[IJK] + 0.5*DTT[IJK]*LL[IJK];

        pgc->start5V(p,F,1);
	}
}





