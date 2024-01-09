/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

void directreini::constraint(lexer *p, fdm* a, ghostcell *pgc, field& b)
{
	
	wallf_update(p,a,pgc);
	
	LOOP
	{
		dirac=0.0;
		
		if(fabs(d0(i,j,k))<epsi)
		dirac = (0.5/epsi)*(1.0 + cos((PI*b(i,j,k))/epsi));
		//dirac = (1.0/epsi)*(1.0 - fabs(d0(i,j,k)));
		
	
		dx = (d0(i+1,j,k) - d0(i-1,j,k))/(2.0*p->DXM); 
		dy = (d0(i,j+1,k) - d0(i,j-1,k))/(2.0*p->DXM); 
		dz = (d0(i,j,k+1) - d0(i,j,k-1))/(2.0*p->DXM); 
		
		dnorm = sqrt(dx*dx + dy*dy + dz*dz);
		
		dnorm = sqrt(dx*dx + dy*dy + dz*dz);
		if(wallf(i,j,k)==1)
		dnorm=1.0;
		
		
		dx = (b(i+1,j,k) - b(i-1,j,k))/(2.0*p->DXM); 
		dy = (b(i,j+1,k) - b(i,j-1,k))/(2.0*p->DXM); 
		dz = (b(i,j,k+1) - b(i,j,k-1))/(2.0*p->DXM); 
		
		
		grad = sqrt(dx*dx + dy*dy + dz*dz);
		if(wallf(i,j,k)==1)
		grad=1.0;
		grad=dnorm;
		
		dval = d0(i,j,k);
		
		sign=dval/sqrt(dval*dval + dnorm*dnorm*p->DXM*p->DXM);
		
		lambda1 =  -dirac * ((b(i,j,k) - d0(i,j,k))/dT) * dV;
		
		lambda2 =  dirac*dirac * (grad) * dV;
		
		
		lambda2 = fabs(lambda2)>1.0e-19?lambda2:1.0e10;
		
		
		b(i,j,k) +=       p->F39*dT*dirac*grad*(lambda1/lambda2);
	}

	
	pgc->start4(p,b,gcval_phi);
}

void directreini::wallf_update(lexer *p, fdm *a, ghostcell *pgc)
{
	int n;
	LOOP
	wallf(i,j,k)=0;
	
	GC4LOOP
	if(p->gcb4[n][4]==21 || p->gcb4[q][4]==22 || p->gcb4[n][4]==5)
	{
	i = p->gcb4[n][0];
	j = p->gcb4[n][1];
	k = p->gcb4[n][2];
	
	wallf(i,j,k)=1;
	}
}
