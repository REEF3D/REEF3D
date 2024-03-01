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

#include"lexer.h"
#include"ghostcell.h"
#include"field.h"
#include"slice.h"

void ghostcell::gcsl_sommerfeld(lexer *p, slice& f, int gcv, int bc, int cs)
{
    double dfx1,dfx4,dfy2,dfy3;
    
    dfx1 = (f(i+1,j)-f(i,j))/p->DXM;
    dfx4 = (f(i,j)-f(i-1,j))/p->DXM;
    dfy2 = (f(i,j)-f(i,j-1))/p->DXM;
    dfy3 = (f(i,j+1)-f(i,j))/p->DXM;

            
	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j) = f(i,j) - p->dt*sqrt(9.81*(p->wd+f(i,j)))*dfx1;

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1) = f(i,j) - p->dt*sqrt(9.81*(p->wd+f(i,j)))*dfy2;

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1) = f(i,j) - p->dt*sqrt(9.81*(p->wd+f(i,j)))*dfy3;

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j) = f(i,j) - p->dt*sqrt(9.81*(p->wd+f(i,j)))*dfx4;

}
