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

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"cpt.h"

void ghostcell::kinematic_fsf(lexer *p,field& f,double dist,int gcv, int bc, int cs)
{
    double wval;
    
	wval = 0.0*(a->eta(i,j) - a->eta_n(i,j))/p->dt
    
         + 0.5*(a->u(i,j,k)+a->u(i-1,j,k))*((a->eta(i+1,j)-a->eta(i-1,j))/(2.0*p->DXP[IP]))
    
         + 0.5*(a->v(i,j,k)+a->v(i,j-1,k))*((a->eta(i,j+1)-a->eta(i,j-1))/(2.0*p->DYP[JP]));

	for(q=0;q<margin;++q)
	f(i,j,k+q+1)= wval;
    
    
    cout<<"KINEMATIC FSF"<<endl;
}

