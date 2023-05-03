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

#include"directreini.h"
#include"lexer.h"
#include"fdm.h"

void directreini::reconstruct(lexer *p,fdm* a, field& b, fieldint& nodeflag, fieldint& vertice)
{
    ccptcount=0;

    for(n=0;n<numtri; ++n)
    confac[n]=-1;
    facount=ccptcount=0;


    for(n=0;n<numtri; ++n)
    {
        if((ls[tri[n][0]] >= -zero && ls[tri[n][1]] < zero)  ||  (ls[tri[n][0]] < zero && ls[tri[n][1]] >= -zero))
        addpoint(p,a,tri[n][0],tri[n][1]);

        if((ls[tri[n][0]] >= -zero && ls[tri[n][2]] < zero)  ||  (ls[tri[n][0]] < zero && ls[tri[n][2]] >= -zero))
        addpoint(p,a,tri[n][0],tri[n][2]);

        if((ls[tri[n][0]] >= -zero && ls[tri[n][3]] < zero)  ||  (ls[tri[n][0]] < zero && ls[tri[n][3]] >= -zero))
        addpoint(p,a,tri[n][0],tri[n][3]);

        if((ls[tri[n][1]] >= -zero && ls[tri[n][2]] < zero)  ||  (ls[tri[n][1]] < zero && ls[tri[n][2]] >= -zero))
        addpoint(p,a,tri[n][1],tri[n][2]);

        if((ls[tri[n][1]] >= -zero && ls[tri[n][3]] < zero)  ||  (ls[tri[n][1]] < zero && ls[tri[n][3]] >= -zero))
        addpoint(p,a,tri[n][1],tri[n][3]);

        if((ls[tri[n][2]] >= -zero && ls[tri[n][3]] < zero)  ||  (ls[tri[n][2]] < zero && ls[tri[n][3]] >= -zero))
        addpoint(p,a,tri[n][2],tri[n][3]);
    }
}

void directreini::addpoint(lexer *p, fdm *a, int q1, int q2)
{
	// p. 917
	
    double dist,xd,dnom;

    //dist = sqrt(pow(pt[q2][0]-pt[q1][0], 2.0) + pow(pt[q2][1]-pt[q1][1], 2.0) + pow(pt[q2][2]-pt[q1][2], 2.0));

    dnom=ls[q2]-ls[q1];
    dnom=fabs(dnom)>1.0e-20?dnom:1.0e-20;

    xd = -(ls[q1]/(dnom));

    ccpt[ccptcount][0] = (pt[q2][0]-pt[q1][0])*xd + pt[q1][0];
    ccpt[ccptcount][1] = (pt[q2][1]-pt[q1][1])*xd + pt[q1][1];
    ccpt[ccptcount][2] = (pt[q2][2]-pt[q1][2])*xd + pt[q1][2];

    if(confac[n]>-1)
    nn=confac[n];

    if(confac[n]==-1)
    {
        confac[n]=facount;
        nn=facount;
        ++facount;
    }

    facet[nn][numfac[n]] = ccptcount;
    ++numfac[n];

    ++ccptcount;
}

