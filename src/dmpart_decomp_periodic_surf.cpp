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

void decomp::periodic_surf(lexer* p, dive* a)
{
    a->periodicXall[0]=a->para1count;
    a->periodicXall[1]=a->para2count;
    a->periodicXall[2]=a->para3count;
    a->periodicXall[3]=a->para4count;
    a->periodicXall[4]=a->para5count;
    a->periodicXall[5]=a->para6count;

    int periodic_count=0;

    // parasurf parallel periodic BC
    if(p->C21==2)
    {
        i=0;
        for(j=0;j<a->knoy;++j)
        for(k=0;k<a->knoz;++k)
        {
            a->para1sf[a->para1count][0]=i;
            a->para1sf[a->para1count][1]=j;
            a->para1sf[a->para1count][2]=k;
            a->para1count++;
            ++periodic_count;
        }

        i=p->knox-1;
        for(j=0;j<a->knoy;++j)
        for(k=0;k<a->knoz;++k)
        {
            a->para4sf[a->para4count][0]=i;
            a->para4sf[a->para4count][1]=j;
            a->para4sf[a->para4count][2]=k;
            a->para4count++;
            ++periodic_count;
        }
    }

    if(p->C22==2)
    {
        j=0;
        for(i=0;i<a->knox;++i)
        for(k=0;k<a->knoz;++k)
        {
            a->para3sf[a->para3count][0]=i;
            a->para3sf[a->para3count][1]=j;
            a->para3sf[a->para3count][2]=k;
            a->para3count++;
            ++periodic_count;
        }

        j=p->knoy-1;
        for(i=0;i<a->knox;++i)
        for(k=0;k<a->knoz;++k)
        {
            a->para2sf[a->para2count][0]=i;
            a->para2sf[a->para2count][1]=j;
            a->para2sf[a->para2count][2]=k;
            a->para2count++;
            ++periodic_count;
        }
    }

    if(p->C23==2)
    {
        k=0;
        for(i=0;i<a->knox;++i)
        for(j=0;j<a->knoy;++j)
        {
            a->para5sf[a->para5count][0]=i;
            a->para5sf[a->para5count][1]=j;
            a->para5sf[a->para5count][2]=k;
            a->para5count++;
            ++periodic_count;
        }

        k=p->knoz-1;
        for(i=0;i<a->knox;++i)
        for(j=0;j<a->knoy;++j)
        {
            a->para6sf[a->para6count][0]=i;
            a->para6sf[a->para6count][1]=j;
            a->para6sf[a->para6count][2]=k;
            a->para6count++;
            ++periodic_count;
        }
    }

    if(p->C21>0||p->C22>0||p->C23>0)
    cout<<"periodic parasf: "<<periodic_count<<endl;
}
}
