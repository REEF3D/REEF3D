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

#include"particle_pls.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void particle_pls::setup(lexer* p, fdm* a, ghostcell* pgc)
{
	LOOP
	{
		posnum(i,j,k)=0.0;
		negnum(i,j,k)=0.0;
	}
	
	pgc->start4(p,posnum,1);
	pgc->start4(p,negnum,1);
	
    dgc_update(p,a,pgc);
    allocate(p,a,pgc);
    reseed(p,a,pgc,1.0);
    setradius(p,a);
    remove(p,a,pgc);
    
	if((p->count%p->F34==0 || p->count==0 )&& p->F34>0)
	{
    print_vtu(p,a,pgc,pos,posflag,posactive,1);
	print_vtu(p,a,pgc,neg,negflag,negactive,2);
	++printcount;
	}
}


