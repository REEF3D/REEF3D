/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::tpflagfield(lexer *p)
{

    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
	p->tpflag[i]=p->flag4[i];

	LOOP
	{
	    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=9;

	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0)
	    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=9;

	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0)
	    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=9;
	}

    LOOP
	{

	    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==9)
	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=11;

	    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==9)
	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=11;

	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==9)
	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==9)
	    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;


	    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==9)
	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==9)
	    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;
	}

    for(n=0;n<p->facetnum;n++)
    {
        i=p->facet[n][0];
        j=p->facet[n][1];
        k=p->facet[n][2];

    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0)
    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=11;

    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=11;

    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0)
    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=11;

    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=11;


    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]<0)
    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;

    if(p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0)
    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=11;

    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]<0)
    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;

    if(p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0)
    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=11;
    }
}

