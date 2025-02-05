/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"strain.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"


double strain::pk(lexer *p, fdm *a,field &eddyv)
{ 
    if(p->j_dir==1)
    {
	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));
    }
    
    if(p->j_dir==0)
    {
	s11 = pudx(p,a);
	s22 = 0.0;
	s33 = pwdz(p,a);
	s12 = 0.0;
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = 0.0;
    }

    val = eddyv(i,j,k)*(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
    
    return val;
    
}

double strain::pk_b(lexer *p, fdm *a,field &eddyv)
{ 
    val =      (1.0/0.85)*(1.0/a->ro(i,j,k))*eddyv(i,j,k)*(
             p->W20*(a->ro(i+1,j,k) - a->ro(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1])
           + p->W21*(a->ro(i,j+1,k) - a->ro(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1])
           + p->W22*(a->ro(i,j,k+1) - a->ro(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));
    
    return val;
}
