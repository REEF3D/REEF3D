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

#include"directreini.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void directreini::correction(lexer *p, fdm* a, ghostcell *pgc, field& b)
{

	double phival,dval,H0,denom;

	for(n=0;n<numvert;++n)
	ls1[n]=ls[n];

		
	for(n=0;n<numvert;++n)
    {
		dV1 = 1.0;
		dV2 = 1.0;
		C1 = 0.0;
		C2 = 0.0;
		mi = 0.0;
		eta = 0.0;
		
		count=0;
		while(dV1>1.0e-10*dV && count<1000)
		{
			
			phival = ls0[n];
			if(phival>epsi)
			H0=1.0;

			if(phival<-epsi)
			H0=0.0;

			if(fabs(phival)<=epsi)
			H0=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
			
			dval = ls[n];
			if(dval>epsi)
			H=1.0;

			if(dval<-epsi)
			H=0.0;

			if(fabs(dval)<=epsi)
			H=0.5*(1.0 + dval/epsi + (1.0/PI)*sin((PI*dval)/epsi));
			
			
			/*
			phival = ls0[n];
			if(phival>=0.0)
			H0=1.0;

			if(phival<0.0)
			H0=0.0;
			
			dval = ls[n];
			if(dval>=0)
			H=1.0;

			if(dval<0)
			H=0.0;
			*/
			
			dV1 = dV*(H0-H);
			
			denom = fabs(ls[n])>1.0e-19?fabs(ls[n]):1.0e20;
			eta = dV1/denom;
			
			ls[n] += eta;
			
			
			++count;
		}
	}
	
	
	for(n=0;n<numvert;++n)
    {
    i=ijk[n][0];
    j=ijk[n][1];
    k=ijk[n][2];
    b(i,j,k)=(1.0-p->F39)*ls1[n] + p->F39*ls[n];
    }

}
