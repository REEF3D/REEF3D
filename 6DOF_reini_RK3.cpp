/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"

void sixdof_f::reini_RK3(lexer* p, fdm* a, ghostcell* pgc, field& b)
{
    n=0;
	LOOP
	{
	dt.V[n] = 0.55*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	++n;
	}

	n=0;
	ALOOP
	{
	f.V[n]=b(i,j,k);
    frk1.V[n]=frk2.V[n]=L.V[n]=0.0;
	++n;
	}
	
    pgc->start4aV(p,f,50);

	reiniter=30;
	
	
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<endl<<"initializing fb..."<<endl<<endl;
	reiniter=25;
	}
    

    for(int q=0;q<reiniter;++q)
    {
	// Step 1
	prdisc->start(p,a,pgc,f,L,5);

	n=0;
    ALOOP
    {
	frk1.V[n] = f.V[n] + dt.V[n]*L.V[n];
    ++n;
    }

	pgc->start4aV(p,frk1,50);

    // Step 2
    prdisc->start(p,a,pgc,frk1,L,5);

	n=0;
    ALOOP
    {
	frk2.V[n]=  0.75*f.V[n] + 0.25*frk1.V[n] + 0.25*dt.V[n]*L.V[n];
    ++n;
    }

	pgc->start4aV(p,frk2,50);

    // Step 3
    prdisc->start(p,a,pgc,frk2,L,5);

	n=0;
    ALOOP
    {
	f.V[n] = (1.0/3.0)*f.V[n] + (2.0/3.0)*frk2.V[n] + (2.0/3.0)*dt.V[n]*L.V[n];
    ++n;
    }

	pgc->start4aV(p,f,50);	
	}
	
	n=0;
	ALOOP
	{
	b(i,j,k)=f.V[n];
	++n;
	}
	
	pgc->start4a(p,b,50);
}



