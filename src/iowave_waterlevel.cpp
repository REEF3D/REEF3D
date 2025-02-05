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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void iowave::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->I230>0)
    ff_waterlevel(p,a,pgc,a->phi);
        
    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void iowave::fsfrkout(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
	for(n=0;n<p->gcout_count;++n)
	{
        i=p->gcout[n][0];
        j=p->gcout[n][1];
        k=p->gcout[n][2];

        f(i+1,j,k)=a->phi(i+1,j,k);
        f(i+2,j,k)=a->phi(i+2,j,k);
        f(i+3,j,k)=a->phi(i+3,j,k);
	}
}

void iowave::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
	for(n=0;n<p->gcin_count;++n)
	{
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
		
        f(i-1,j,k)=a->phi(i-1,j,k);
        f(i-2,j,k)=a->phi(i-2,j,k);
        f(i-3,j,k)=a->phi(i-3,j,k);
    }
}

