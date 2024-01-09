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

void ghostcell::flagfield(lexer *p)
{
    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
    p->flag[i]=1;

        
    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
    {
    if(p->flag4[i]==1)
    p->flag4[i]=10;

    if(p->flag4[i]==-1)
    p->flag4[i]=OBJ;
    }
    
    flagx(p,p->flag4);
    
	if(p->Y60==1)
    LOOP
    {   
        if(p->i_dir==1)
        if(p->flag4[Im1JK]<0
        && p->flag4[Ip1JK]<0)
        p->flag4[IJK]=OBJ;
        
        if(p->j_dir==1)
        if(p->flag4[IJm1K]<0
        && p->flag4[IJp1K]<0)
        p->flag4[IJK]=OBJ;
        
        if(p->k_dir==1)
        if(p->flag4[IJKm1]<0
        && p->flag4[IJKp1]<0)
        p->flag4[IJK]=OBJ;
    }
    


    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
	{
	p->flag1[i]=p->flag4[i];
	p->flag2[i]=p->flag4[i];
	p->flag3[i]=p->flag4[i];
	}

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==4 && (p->periodic1!=1 || i+p->origin_i<p->gknox-1))
        p->flag1[IJK]=OBJ;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==2 && (p->periodic2!=1 || j+p->origin_j<p->gknoy-1))
        p->flag2[IJK]=OBJ;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==6 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
        p->flag3[IJK]=OBJ;
    }
	
}

void ghostcell::flagfield_topo(lexer *p)
{

    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
    {
    if(p->flag4[i]==1)
    p->flag4[i]=10;

    if(p->flag4[i]==-1)
    p->flag4[i]=OBJ;
    }
    
    flagx(p,p->flag4);
    
	if(p->Y60==1)
    LOOP
    {
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0
        && p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=OBJ;

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0
        && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]<0)
        p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=OBJ;

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0
        && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]<0)
        p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=OBJ;
    }


    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
	{
	p->flag1[i]=p->flag4[i];
	p->flag2[i]=p->flag4[i];
	p->flag3[i]=p->flag4[i];
	}

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==4)
        p->flag1[IJK]=OBJ;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==2)
        p->flag2[IJK]=OBJ;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==6)
        p->flag3[IJK]=OBJ;
    }
	
}

