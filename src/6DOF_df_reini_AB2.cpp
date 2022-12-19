/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"

void sixdof_df_object::reini_AB2(lexer* p, fdm* a, ghostcell* pgc, field& b)
{	
	n=0;
	ALOOP
	{
	f.V[n]=b(i,j,k);
	++n;
	}
    
    pgc->start4avec(p,f,50);
	
    n=0;
	LOOP
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
		prdisc->start(p,a,pgc,f,L,5);

		if(q==0)
		NLOOP4A
		frk1.V[n]=L.V[n];


		NLOOP4A
		{
		f.V[n] += dt.V[n]*0.5*(3.0*L.V[n] - frk1.V[n]);

		frk1.V[n]=L.V[n];
		}

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





