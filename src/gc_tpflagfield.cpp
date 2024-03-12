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

void ghostcell::tpflagfield(lexer *p)
{    
    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
    p->tpflag[i]=1;
    
    LOOP
	p->tpflag[IJK]=p->flag4[IJK];

	LOOP
	{
	    if(p->tpflag[Im1JK]<0)
	    p->tpflag[Im1JK]=9;

	    if(p->tpflag[IJm1K]<0)
	    p->tpflag[IJm1K]=9;

	    if(p->tpflag[IJKm1]<0)
	    p->tpflag[IJKm1]=9;
	}

    LOOP
	{

	    if(p->tpflag[Im1JK]==9)
	    if(p->tpflag[IJm1K]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=11;

	    if(p->tpflag[Im1JK]==9)
	    if(p->tpflag[IJKm1]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=11;

	    if(p->tpflag[IJm1K]==9)
	    if(p->tpflag[IJKm1]==9)
	    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;


	    if(p->tpflag[Im1JK]==9)
	    if(p->tpflag[IJm1K]==9)
	    if(p->tpflag[IJKm1]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;
	}
}

void ghostcell::tpflagfield_sigma(lexer *p)
{    
    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
    p->tpflag[i]=1;
    /*
    LOOP
	p->tpflag[IJK]=p->flag4[IJK];

	LOOP
	{
	    if(p->tpflag[Im1JK]<0)
	    p->tpflag[Im1JK]=9;

	    if(p->tpflag[IJm1K]<0)
	    p->tpflag[IJm1K]=9;

	    if(p->tpflag[IJKm1]<0)
	    p->tpflag[IJKm1]=9;
	}

    LOOP
	{

	    if(p->tpflag[Im1JK]==9)
	    if(p->tpflag[IJm1K]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]=11;

	    if(p->tpflag[Im1JK]==9)
	    if(p->tpflag[IJKm1]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]=11;

	    if(p->tpflag[IJm1K]==9)
	    if(p->tpflag[IJKm1]==9)
	    p->tpflag[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;


	    if(p->tpflag[Im1JK]==9)
	    if(p->tpflag[IJm1K]==9)
	    if(p->tpflag[IJKm1]==9)
	    p->tpflag[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1]=11;
	}*/
}
