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

#include"sflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sflow_f::maxcoor(lexer *p, fdm2D* b, ghostcell* pgc)
{
	p->maxlength=-1.0e9;
	p->xcoormax=-1.0e9;
	p->xcoormin=1.0e9;
	p->ycoormax=-1.0e9;
	p->ycoormin=1.0e9;
	p->zcoormax=-1.0e9;
	p->zcoormin=1.0e9;

    SLICELOOP4
    {
        p->xcoormax = MAX(p->xcoormax,p->XN[IP1]);
        p->xcoormin = MIN(p->xcoormin,p->XN[IP]);
        p->ycoormax = MAX(p->ycoormax,p->YN[JP1]);
        p->ycoormin = MIN(p->ycoormin,p->YN[JP]);
     }

     p->maxlength=MAX(p->maxlength,p->xcoormax-p->xcoormin);
     p->maxlength=MAX(p->maxlength,p->ycoormax-p->ycoormin);

     p->maxlength=pgc->globalmax(p->maxlength);
	 
	 p->xcoormax=pgc->globalmax(p->xcoormax);
	 p->ycoormax=pgc->globalmax(p->ycoormax);
	 
	 p->xcoormin=pgc->globalmin(p->xcoormin);
	 p->ycoormin=pgc->globalmin(p->ycoormin);
	 
}
