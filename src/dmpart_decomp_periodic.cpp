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

void decomp::periodic_ini(lexer* p, dive* a)
{
    // serial or parallel periodic BC
    if(p->C21>=1)
    {
        if(a->mx==1)
        p->C21=1;
        else if(a->mx>1)
        p->C21=2;
    }

    if(p->C22>=1)
    {
        if(a->my==1)
        p->C22=1;
        else if(a->my>1)
        p->C22=2;
    }

    if(p->C23>=1)
    {
        if(a->mz==1)
        p->C23=1;
        else if(a->mz>1)
        p->C23=2;
    }
}

void decomp::periodic_nb(lexer* p, dive* a)
{
    // x-dir parallel periodic BC
    if(p->C21==2)
    {
        count=0;
        NLOOP
        {
            count++;

            if(aa==1)
            a->nbpara1[count]=a->sgfield[a->mx][bb][cc]-1;

            if(aa==a->mx)
            a->nbpara4[count]=a->sgfield[1][bb][cc]-1;
        }
    }

    // y-dir parallel periodic BC
    if(p->C22==2)
    {
        count=0;
        NLOOP
        {
            count++;

            if(bb==1)
            a->nbpara3[count]=a->sgfield[aa][a->my][cc]-1;

            if(bb==a->my)
            a->nbpara2[count]=a->sgfield[aa][1][cc]-1;
        }
    }

    // z-dir parallel periodic BC
    if(p->C23==2)
    {
        count=0;
        NLOOP
        {
            count++;

            if(cc==1)
            a->nbpara5[count]=a->sgfield[aa][bb][a->mz]-1;

            if(cc==a->mz)
            a->nbpara6[count]=a->sgfield[aa][bb][1]-1;
        }
    }
}

void decomp::periodic_count(lexer* p, dive* a)
{
    for(i=0;i<a->periodicXall[0];i++)
    {
        n=a->subgrid(a->para1sf[i][0],a->para1sf[i][1],a->para1sf[i][2]);
        a->periodicX[n][0]++;
    }

    for(i=0;i<a->periodicXall[1];i++)
    {
        n=a->subgrid(a->para2sf[i][0],a->para2sf[i][1],a->para2sf[i][2]);
        a->periodicX[n][1]++;
    }

    for(i=0;i<a->periodicXall[2];i++)
    {
        n=a->subgrid(a->para3sf[i][0],a->para3sf[i][1],a->para3sf[i][2]);
        a->periodicX[n][2]++;
    }

    for(i=0;i<a->periodicXall[3];i++)
    {
        n=a->subgrid(a->para4sf[i][0],a->para4sf[i][1],a->para4sf[i][2]);
        a->periodicX[n][3]++;
    }

    for(i=0;i<a->periodicXall[4];i++)
    {
        n=a->subgrid(a->para5sf[i][0],a->para5sf[i][1],a->para5sf[i][2]);
        a->periodicX[n][4]++;
    }

    for(i=0;i<a->periodicXall[5];i++)
    {
        n=a->subgrid(a->para6sf[i][0],a->para6sf[i][1],a->para6sf[i][2]);
        a->periodicX[n][5]++;
    }
}
}
