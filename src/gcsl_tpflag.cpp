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

void ghostcell::gcsl_tpflag(lexer *p)
{
	for(i=0;i<imax*jmax;++i)
	p->tpflagslice[i]=-1;
	
	
	SLICEBASELOOP
	p->tpflagslice[IJ]=1;

	SLICEBASELOOP
	{
	    if(p->tpflagslice[Im1J]<=0)
	    p->tpflagslice[Im1J]=9;

	    if(p->tpflagslice[IJm1]<=0)
	    p->tpflagslice[IJm1]=9;

	}

    SLICEBASELOOP
	{
	    if(p->tpflagslice[Im1J]==9)
	    if(p->tpflagslice[IJm1]==9)
	    p->tpflagslice[Im1Jm1]=11;
	}
	
	/*
    for(i=0;i<p->imax*p->jmax; ++i)
	p->tpflagslice[i]=p->flagslice4[i];

	SLICELOOP4
	{
	    if(p->tpflagslice[Im1J]<=0)
	    p->tpflagslice[Im1J]=9;

	    if(p->tpflagslice[IJm1]<=0)
	    p->tpflagslice[IJm1]=9;

	}

    SLICELOOP4
	{
	    if(p->tpflagslice[Im1J]==9)
	    if(p->tpflagslice[IJm1]==9)
	    p->tpflagslice[Im1Jm1]=11;
	}*/
}

