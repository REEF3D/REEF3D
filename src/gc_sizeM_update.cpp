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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::sizeM_update(lexer* p, fdm* a)
{
	count=0;
	ULOOP
	++count;
	
	p->sizeM1[0]=0;
	p->sizeM1[1]=count;
	
	count=0;
	VLOOP
	++count;
	
	p->sizeM2[0]=0;
	p->sizeM2[1]=count;
	
	count=0;
	WLOOP
	++count;
	
	p->sizeM3[0]=0;
	p->sizeM3[1]=count;
	
	count=0;
	LOOP
	++count;
	
	p->sizeM4[0]=0;
	p->sizeM4[1]=count;
    
    count=0;
	ALOOP
	++count;
	
	p->sizeM4a[0]=0;
	p->sizeM4a[1]=count;
    
    
    count=0;
	BASELOOP
	++count;
	
	p->sizeM6[0]=0;
	p->sizeM6[1]=count;
    
    p->sizeM9[0]=0;
	p->sizeM9[1]=count;
	
}
