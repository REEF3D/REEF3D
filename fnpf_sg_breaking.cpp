/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_sg_fsfbc.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_sg_fsfbc::breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, double alpha)
{
    if(p->A346==1)
    SLICELOOP4
    {
            diss(i,j)=0.0;
            
            if( (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A347*sqrt(9.81*c->WL(i,j)))
            diss(i,j)=1.86;
    }
    
    
    
    
    
}