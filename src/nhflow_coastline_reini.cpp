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

#include"nhflow_coastline.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

void nhflow_coastline::reini(lexer *p, ghostcell *pgc, slice &f)
{
	if(p->count==0)
	{
	reiniter=2*int(p->maxlength/(1.0*p->DXM));
    
    if(p->mpirank==0)
	cout<<"initializing coastline... "<<endl<<endl;
	}

	if(p->count>0)
	step(p);

    for(int q=0;q<reiniter;++q)
    {
	// Step 1
    disc(p,pgc,f);
    
	SLICELOOP4
	frk1(i,j) = f(i,j) + dt(i,j)*L(i,j);

	pgc->gcsl_start4(p,frk1,50);
    

    // Step 2
    disc(p,pgc,frk1);
    
	SLICELOOP4
	frk2(i,j)=  0.75*f(i,j) + 0.25*frk1(i,j) + 0.25*dt(i,j)*L(i,j);

	pgc->gcsl_start4(p,frk2,50);


    // Step 3
    disc(p,pgc,frk2);
    
	SLICELOOP4
	f(i,j) = (1.0/3.0)*f(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*dt(i,j)*L(i,j);

	pgc->gcsl_start4(p,f,50);
	}
    
    pgc->gcsl_start4(p,f,50);
}


void nhflow_coastline::step(lexer* p)
{
	reiniter=p->S37;
}

void nhflow_coastline::time_preproc(lexer* p)
{	
    n=0;
	SLICELOOP4
	{
	dt(i,j) = p->F43*MIN(p->DXP[IP],p->DYP[JP]);
	++n;
	}
}

