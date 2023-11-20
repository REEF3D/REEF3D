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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"

void sixdof_df_object::reini_RK2(lexer* p, fdm* a, ghostcell* pgc, field& b)
{	
	n=0;
	ALOOP
	{
	f.V[n]=b(i,j,k);
	++n;
	}
    
    pgc->start4avec(p,f,50);
	
    n=0;
	ALOOP
	{
	dt.V[n] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	++n;
	}
	
	reiniter=10;
	
	
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<endl<<"initializing fb..."<<endl<<endl;
	reiniter=10;
	}

    for(int q=0;q<reiniter;++q)
	{
        // Step 1
		prdisc->start(p,a,pgc,f,L,5);

		NLOOP4A
		frk1.V[n] = f.V[n] + dt.V[n]*L.V[n];

         pgc->start4avec(p,frk1,50);
        
        
        // Step 2
		prdisc->start(p,a,pgc,frk1,L,5);

		NLOOP4A
		f.V[n] = 0.5*f.V[n] + 0.5*frk1.V[n] + 0.5*dt.V[n]*L.V[n];

        pgc->start4avec(p,f,50);
	}
		
	n=0;
	ALOOP
	{
	b(i,j,k)=f.V[n];
	++n;
	}
	
	pgc->start4a(p,b,50);
}





