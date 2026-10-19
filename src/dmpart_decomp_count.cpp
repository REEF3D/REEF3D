/*--------------------------------------------------------------------
DIVEMesh
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D (port of the DIVEMesh decomposition, dmpart).

DIVEMesh is free software; you can redistribute it and/or modify it
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

#include "dmpart_decomp.h"

namespace dmpart
{

void decomp::surfcount(lexer* p,dive* a)
{
    for(i=0;i<a->para1count;i++)
    {
        n=a->subgrid(a->para1sf[i][0],a->para1sf[i][1],a->para1sf[i][2]);
        a->para1[n]++;
    }

    for(i=0;i<a->para2count;i++)
    {
        n=a->subgrid(a->para2sf[i][0],a->para2sf[i][1],a->para2sf[i][2]);
        a->para2[n]++;
    }

    for(i=0;i<a->para3count;i++)
    {
        n=a->subgrid(a->para3sf[i][0],a->para3sf[i][1],a->para3sf[i][2]);
        a->para3[n]++;
    }

    for(i=0;i<a->para4count;i++)
    {
        n=a->subgrid(a->para4sf[i][0],a->para4sf[i][1],a->para4sf[i][2]);
        a->para4[n]++;
    }

    for(i=0;i<a->para5count;i++)
    {
        n=a->subgrid(a->para5sf[i][0],a->para5sf[i][1],a->para5sf[i][2]);
        a->para5[n]++;
    }

    for(i=0;i<a->para6count;i++)
    {
        n=a->subgrid(a->para6sf[i][0],a->para6sf[i][1],a->para6sf[i][2]);
        a->para6[n]++;
    }
}

void decomp::slicesurfcount(lexer* p,dive* a)
{
    for(i=0;i<a->paraslice1count;i++)
    {
        n=a->subslice(a->paraslice1sf[i][0],a->paraslice1sf[i][1]);
        a->paraslice1[n]++;
    }

    for(i=0;i<a->paraslice2count;i++)
    {
        n=a->subslice(a->paraslice2sf[i][0],a->paraslice2sf[i][1]);
        a->paraslice2[n]++;
    }

    for(i=0;i<a->paraslice3count;i++)
    {
        n=a->subslice(a->paraslice3sf[i][0],a->paraslice3sf[i][1]);
        a->paraslice3[n]++;
    }

    for(i=0;i<a->paraslice4count;i++)
    {
        n=a->subslice(a->paraslice4sf[i][0],a->paraslice4sf[i][1]);
        a->paraslice4[n]++;
    }
}

void decomp::slicecornercount(lexer* p,dive* a)
{
    for(i=0;i<a->paracoslice1count;i++)
    {
        n=a->paracoslice1sf[i][2];
        a->paracoslice1[n]++;
    }

    for(i=0;i<a->paracoslice2count;i++)
    {
        n=a->paracoslice2sf[i][2];
        a->paracoslice2[n]++;
    }

    for(i=0;i<a->paracoslice3count;i++)
    {
        n=a->paracoslice3sf[i][2];
        a->paracoslice3[n]++;
    }

    for(i=0;i<a->paracoslice4count;i++)
    {
        n=a->paracoslice4sf[i][2];
        a->paracoslice4[n]++;
    }
}

void decomp::cornercount(lexer* p,dive* a)
{

    for(i=0;i<a->paraco1count;i++)
    {
        n=a->para1co[i][3];
        a->paraco1[n]++;
    }

    for(i=0;i<a->paraco2count;i++)
    {
        n=a->para2co[i][3];
        a->paraco2[n]++;
    }

    for(i=0;i<a->paraco3count;i++)
    {
        n=a->para3co[i][3];
        a->paraco3[n]++;
    }

    for(i=0;i<a->paraco4count;i++)
    {
        n=a->para4co[i][3];
        a->paraco4[n]++;
    }

    for(i=0;i<a->paraco5count;i++)
    {
        n=a->para5co[i][3];
        a->paraco5[n]++;
    }

    for(i=0;i<a->paraco6count;i++)
    {
        n=a->para6co[i][3];
        a->paraco6[n]++;
    }
}
}
