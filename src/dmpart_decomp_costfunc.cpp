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

void decomp::nodecalc(lexer *p, dive *a)
{
    numfactor=0;

    for(i=1;i<=p->M10;i++)
    for(j=1;j<=p->M10;j++)
    for(k=1;k<=p->M10;k++)
    if(i*j*k==p->M10)
    numfactor++;


    p->Iarray(a->mpx,numfactor+1);
    p->Iarray(a->mpy,numfactor+1);
    p->Iarray(a->mpz,numfactor+1);

    numfactor=0;

    for(i=1;i<=p->M10;i++)
    for(j=1;j<=p->M10;j++)
    for(k=1;k<=p->M10;k++)
    if(i*j*k==p->M10)
    {
        a->mpx[numfactor]=i;
        a->mpy[numfactor]=j;
        a->mpz[numfactor]=k;
        numfactor++;
    }
}

void decomp::costfunc(lexer* p, dive* a)
{
    int n;
    costmin=1.0e18;

    a->mx=a->my=a->mz=0;

    for(n=0;n<numfactor;n++)
    {
        if(a->mpx[n]>a->knox || a->mpy[n]>a->knoy || a->mpz[n]>a->knoz)
        continue;

        cost = double(a->knox)/double(a->mpx[n])*double(a->knoy)/double(a->mpy[n])
             + double(a->knox)/double(a->mpx[n])*double(a->knoz)/double(a->mpz[n])
             + double(a->knoy)/double(a->mpy[n])*double(a->knoz)/double(a->mpz[n]);

        if(costmin>=cost && p->M11==1 && p->M12==1 && p->M13==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }

        if(costmin>=cost && p->M11==0 && p->M12==1 && p->M13==1 && a->mpx[n]==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }

        if(costmin>=cost && p->M11==1 && p->M12==0 && p->M13==1 && a->mpy[n]==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }

        if(costmin>=cost && p->M11==1 && p->M12==1 && p->M13==0 && a->mpz[n]==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }

        if(costmin>=cost && p->M11==1 && p->M12==0 && p->M13==0 && a->mpy[n]==1 && a->mpz[n]==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }

        if(costmin>=cost && p->M11==0 && p->M12==1 && p->M13==0 && a->mpx[n]==1 && a->mpz[n]==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }

        if(costmin>=cost && p->M11==0 && p->M12==0 && p->M13==1 && a->mpx[n]==1 && a->mpy[n]==1)
        {
            costmin=cost;
            a->mx=a->mpx[n];
            a->my=a->mpy[n];
            a->mz=a->mpz[n];
        }
    }

    if(a->mx*a->my*a->mz!=p->M10)
    {
        cout<<endl<<"!!! no domain decomposition for M 10 "<<p->M10<<" fits the grid ("
            <<a->knox<<" x "<<a->knoy<<" x "<<a->knoz<<" cells) with the directions enabled by M 11/12/13 !!!"<<endl<<endl;
        exit(1);
    }

    ddout<<"partition:  "<<a->mx<<' '<<a->my<<' '<<a->mz<<endl;
}
}
