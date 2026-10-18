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

void decomp::paracosurface(lexer* p, dive* a)
{
    int r;
    int paracell[3][3][3];

    // CORNERS
    a->paraco1count=0;
    a->paraco2count=0;
    a->paraco3count=0;
    a->paraco4count=0;
    a->paraco5count=0;
    a->paraco6count=0;

    LOOP
    {

        for(n=0;n<3;++n)
        for(q=0;q<3;++q)
        for(r=0;r<3;++r)
        paracell[n][q][r]=0;

        //  1
        if(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i-1,j,k))))
        {
            paracell[0][0][1]+=1;
        }

        if(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i-1,j,k))))
        {
            paracell[0][2][1]+=1;
        }

        if(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i-1,j,k))))
        {
            paracell[0][1][0]+=1;
        }

        if(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i-1,j,k))))
        {
            paracell[0][1][2]+=1;
        }

        if(p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[0][0][1] += 1;
        }

        if(p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[0][2][1] += 1;
        }

        if(p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) )))
        {
            paracell[0][1][0] += 1;
        }

        if(p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) )))
        {
            paracell[0][1][2] += 1;
        }


        //  4
        if(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i+1,j,k))))
        paracell[2][0][1]+=1;

        if(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i+1,j,k))))
        paracell[2][2][1]+=1;

        if(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i+1,j,k))))
        paracell[2][1][0]+=1;

        if(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i+1,j,k))))
        paracell[2][1][2]+=1;

        if(p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[2][0][1] += 1;
        }


        if(p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[2][2][1] += 1;
        }

        if(p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) )))
        {
            paracell[2][1][0] += 1;
        }

        if(p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) )))
        {
            paracell[2][1][2] += 1;
        }

        //  3
        if(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i,j-1,k))))
        paracell[0][0][1]+=1;

        if(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i,j-1,k))))
        paracell[2][0][1]+=1;

        if(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j-1,k))))
        paracell[1][0][0]+=1;

        if(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j-1,k))))
        paracell[1][0][2]+=1;

        if(p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[0][0][1] += 1;
        }

        if(p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[2][0][1] += 1;
        }

        if(p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) )))
        {
            paracell[1][0][0] += 1;
        }

        if(p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) )))
        {
            paracell[1][0][2] += 1;
        }

        //  2
        if(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i,j+1,k))))
        paracell[0][2][1]+=1;

        if(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i,j+1,k))))
        paracell[2][2][1]+=1;

        if(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j+1,k))))
        paracell[1][2][0]+=1;

        if(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j+1,k))))
        paracell[1][2][2]+=1;

        if(p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[0][2][1] += 1;
        }

        if(p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[2][2][1] += 1;
        }

        if(p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) )))
        {
            paracell[1][2][0] += 1;
        }

        if(p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) )))
        {
            paracell[1][2][2] += 1;
        }

        //  5
        if(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i,j,k-1))))
        paracell[0][1][0]+=1;

        if(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i,j,k-1))))
        paracell[2][1][0]+=1;

        if(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j,k-1))))
        paracell[1][0][0]+=1;

        if(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j,k-1))))
        paracell[1][2][0]+=1;

        if(p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[0][1][0] += 1;
        }

        if(p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[2][1][0] += 1;
        }

        if(p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[1][0][0] += 1;
        }

        if(p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[1][2][0] += 1;
        }

        //  6
        if(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i,j,k+1))))
        paracell[0][1][2]+=1;

        if(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i,j,k+1))))
        paracell[2][1][2]+=1;

        if(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j,k+1))))
        paracell[1][0][2]+=1;

        if(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j,k+1))))
        paracell[1][2][2]+=1;

        if(p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[0][1][2] += 1;
        }

        if(p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) )))
        {
            paracell[2][1][2] += 1;
        }

        if(p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[1][0][2] += 1;
        }

        if(p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) )))
        {
            paracell[1][2][2] += 1;
        }

        //-----------------------------------------------------------------------------------------------
        //  1
        if((a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i-1,j,k))))
        || (p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j-1;
            a->para1co[a->paraco1count][2]=k;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=3;
            a->para1co[a->paraco1count][6]=paracell[0][0][1];
            ++a->paraco1count;
        }

        if((a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i-1,j,k))))
        || (p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j+1;
            a->para1co[a->paraco1count][2]=k;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=2;
            a->para1co[a->paraco1count][6]=paracell[0][2][1];
            a->paraco1count++;
        }

        if((a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i-1,j,k))))
        || (p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) ))))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j;
            a->para1co[a->paraco1count][2]=k-1;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=5;
            a->para1co[a->paraco1count][6]=paracell[0][1][0];
            a->paraco1count++;
        }

        if((a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i-1,j,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i-1,j,k))))
        || (p->C21 == 2
        && a->flag(i-1,j,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  ))))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j;
            a->para1co[a->paraco1count][2]=k+1;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=6;
            a->para1co[a->paraco1count][6]=paracell[0][1][2];
            a->paraco1count++;
        }

        //  4
        if((a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i+1,j,k))))
        || (p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j-1;
            a->para4co[a->paraco4count][2]=k;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=3;
            a->para4co[a->paraco4count][6]=paracell[2][0][1];
            a->paraco4count++;
        }


        if((a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i+1,j,k))))
        || (p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j+1;
            a->para4co[a->paraco4count][2]=k;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=2;
            a->para4co[a->paraco4count][6]=paracell[2][2][1];
            a->paraco4count++;
        }

        if((a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i+1,j,k))))
        || (p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) ))))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j;
            a->para4co[a->paraco4count][2]=k-1;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=5;
            a->para4co[a->paraco4count][6]=paracell[2][1][0];
            a->paraco4count++;
        }

        if((a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i+1,j,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i+1,j,k))))
        || (p->C21 == 2
        && a->flag(i+1,j,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) ))))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j;
            a->para4co[a->paraco4count][2]=k+1;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=6;
            a->para4co[a->paraco4count][6]=paracell[2][1][2];
            a->paraco4count++;
        }

        //  3
        if((a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i,j-1,k))))
        || (p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j-1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para3co[a->paraco3count][0]=i-1;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=1;
            a->para3co[a->paraco3count][6]=paracell[0][0][1];
            a->paraco3count++;
        }

        if((a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i,j-1,k))))
        || (p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j-1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para3co[a->paraco3count][0]=i+1;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=4;
            a->para3co[a->paraco3count][6]=paracell[2][0][1];
            a->paraco3count++;
        }

        if((a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j-1,k))))
        || (p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) ))))
        {
            a->para3co[a->paraco3count][0]=i;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k-1;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=5;
            a->para3co[a->paraco3count][6]=paracell[1][0][0];
            a->paraco3count++;
        }

        if((a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j-1,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j-1,k))))
        || (p->C22 == 2
        && a->flag(i,j-1,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) ))))
        {
            a->para3co[a->paraco3count][0]=i;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k+1;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=6;
            a->para3co[a->paraco3count][6]=paracell[1][0][2];
            a->paraco3count++;
        }

        //  2
        if((a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i,j+1,k))))
        || (p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j+1,k)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para2co[a->paraco2count][0]=i-1;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=1;
            a->para2co[a->paraco2count][6]=paracell[0][2][1];
            a->paraco2count++;
        }

        if((a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i,j+1,k))))
        || (p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j+1,k)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para2co[a->paraco2count][0]=i+1;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=4;
            a->para2co[a->paraco2count][6]=paracell[2][2][1];
            a->paraco2count++;
        }

        if((a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j+1,k))))
        || (p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i,j,k-1)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j,k-1)!=a->subgrid(i,j,k) ))))
        {
            a->para2co[a->paraco2count][0]=i;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k-1;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=5;
            a->para2co[a->paraco2count][6]=paracell[1][2][0];
            a->paraco2count++;
        }

        if((a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j+1,k)>0
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j+1,k))))
        || (p->C22 == 2
        && a->flag(i,j+1,k) == -21
        && ((a->flag(i,j,k+1)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j,k+1)!=a->subgrid(i,j,k) ))))
        {
            a->para2co[a->paraco2count][0]=i;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k+1;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=6;
            a->para2co[a->paraco2count][6]=paracell[1][2][2];
            a->paraco2count++;
        }

        //  5
        if((a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i,j,k-1))))
        || (p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k-1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para5co[a->paraco5count][0]=i-1;
            a->para5co[a->paraco5count][1]=j;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=1;
            a->para5co[a->paraco5count][6]=paracell[0][1][0];
            a->paraco5count++;
        }

        if((a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i,j,k-1))))
        || (p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k-1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para5co[a->paraco5count][0]=i+1;
            a->para5co[a->paraco5count][1]=j;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=4;
            a->para5co[a->paraco5count][6]=paracell[2][1][0];
            a->paraco5count++;
        }

        if((a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j,k-1))))
        || (p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k-1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para5co[a->paraco5count][0]=i;
            a->para5co[a->paraco5count][1]=j-1;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=3;
            a->para5co[a->paraco5count][6]=paracell[1][0][0];
            a->paraco5count++;
        }

        if((a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k-1)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j,k-1))))
        || (p->C23 == 2
        && a->flag(i,j,k-1) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k-1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para5co[a->paraco5count][0]=i;
            a->para5co[a->paraco5count][1]=j+1;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=2;
            a->para5co[a->paraco5count][6]=paracell[1][2][0];
            a->paraco5count++;
        }

        //  6
        if((a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i,j,k+1))))
        || (p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i-1,j,k)<0  && a->flag(i-1,j,k+1)<0)
            ||(a->subgrid(i-1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para6co[a->paraco6count][0]=i-1;
            a->para6co[a->paraco6count][1]=j;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=1;
            a->para6co[a->paraco6count][6]=paracell[0][1][2];
            a->paraco6count++;
        }

        if((a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i,j,k+1))))
        || (p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i+1,j,k)<0  && a->flag(i+1,j,k+1)<0)
            ||(a->subgrid(i+1,j,k)!=a->subgrid(i,j,k) ))))
        {
            a->para6co[a->paraco6count][0]=i+1;
            a->para6co[a->paraco6count][1]=j;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=4;
            a->para6co[a->paraco6count][6]=paracell[2][1][2];
            a->paraco6count++;
        }

        if((a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j,k+1))))
        || (p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i,j-1,k)<0  && a->flag(i,j-1,k+1)<0)
            ||(a->subgrid(i,j-1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para6co[a->paraco6count][0]=i;
            a->para6co[a->paraco6count][1]=j-1;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=3;
            a->para6co[a->paraco6count][6]=paracell[1][0][2];
            a->paraco6count++;
        }

        if((a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)
        && a->flag(i,j,k)>0
        && a->flag(i,j,k+1)>0
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j,k+1))))
        || (p->C23 == 2
        && a->flag(i,j,k+1) == -21
        && ((a->flag(i,j+1,k)<0  && a->flag(i,j+1,k+1)<0)
            ||(a->subgrid(i,j+1,k)!=a->subgrid(i,j,k) ))))
        {
            a->para6co[a->paraco6count][0]=i;
            a->para6co[a->paraco6count][1]=j+1;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=2;
            a->para6co[a->paraco6count][6]=paracell[1][2][2];
            a->paraco6count++;
        }
    }
}

void decomp::cornersurface(lexer* p, dive* a)
{
    bool condition = 0;
    bool condperiod = 0;
    bool isupper, isright, islower, isleft;
    bool isupperperiod, isrightperiod, islowerperiod, isleftperiod;
    int r;
    int paracell[3][3][3];

    LOOP
    {
        for(n=0;n<3;++n)
        for(q=0;q<3;++q)
        for(r=0;r<3;++r)
        paracell[n][q][r]=0;

        //  1
        condition =  a->flag(i,j,k)>0 && a->flag(i-1,j,k)>0 && a->subgrid(i,j,k) != a->subgrid(i-1,j,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i-1,j,k);
        isright   = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i-1,j,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i-1,j,k);
        isleft    = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i-1,j,k);

        condperiod       = ( p->C21 == 2) && (a->flag(i-1,j,k) == -21);
        isupperperiod    = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isrightperiod    = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod    = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleftperiod     = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        paracell[0][2][0]+=1;

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        paracell[0][2][2]+=1;

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        paracell[0][0][0]+=1;

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        paracell[0][0][2]+=1;

        //  4
        condition = a->flag(i,j,k)>0 && a->flag(i+1,j,k)>0 && a->subgrid(i,j,k) != a->subgrid(i+1,j,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i+1,j,k);
        isright   = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i+1,j,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i+1,j,k);
        isleft    = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i+1,j,k);

        condperiod = (p->C21 == 2) && (a->flag(i+1,j,k) == -21);
        isupperperiod   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isrightperiod   = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleftperiod    = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        paracell[2][2][0]+=1;

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        paracell[2][2][2]+=1;

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        paracell[2][0][0]+=1;

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        paracell[2][0][2]+=1;

        //  3
        condition = a->flag(i,j,k)>0 && a->flag(i,j-1,k)>0 && a->subgrid(i,j,k) != a->subgrid(i,j-1,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j-1,k);
        isright   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i,j-1,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j-1,k);
        isleft    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i,j-1,k);

        condperiod = (p->C22 == 2) && (a->flag(i,j-1,k) == -21);
        isupperperiod   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isrightperiod   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);
        islowerperiod   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleftperiod    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        paracell[0][0][0]+=1;

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        paracell[0][0][2]+=1;

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        paracell[2][0][0]+=1;

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        paracell[2][0][2]+=1;

        //  2
        condition = a->flag(i,j,k)>0 && a->flag(i,j+1,k)>0 && a->subgrid(i,j,k) != a->subgrid(i,j+1,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j+1,k);
        isright   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i,j+1,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j+1,k);
        isleft    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i,j+1,k);

        condperiod = (p->C22 == 2) && (a->flag(i,j+1,k) == -21);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isright   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleft    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        paracell[0][2][0]+=1;

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        paracell[0][2][2]+=1;

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        paracell[2][2][0]+=1;

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        paracell[2][2][2]+=1;

        //  5
        condition = a->flag(i,j,k)>0 && a->flag(i,j,k-1)>0 && a->subgrid(i,j,k) != a->subgrid(i,j,k-1);
        isupper = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j,k-1);
        islower  = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j,k-1);
        isleft  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i,j,k-1);
        isright = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i,j,k-1);

        condperiod = (p->C23 == 2) && (a->flag(i,j,k-1) == -21);
        isupperperiod = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);
        isleftperiod  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);
        isrightperiod = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);

        if((condition && islower && isleft)||(condperiod && islowerperiod && isleftperiod))
        paracell[0][0][0]+=1;

        if((condition && isupper && isleft)||(condperiod && isupperperiod && isleftperiod))
        paracell[0][2][0]+=1;

        if((condition && islower && isright)||(condperiod && islowerperiod && isrightperiod))
        paracell[2][0][0]+=1;

        if((condition && isupper && isright)||(condperiod && isupperperiod && isrightperiod))
        paracell[2][2][0]+=1;

        //  6
        condition = a->flag(i,j,k)>0 && a->flag(i,j,k+1)>0 && a->subgrid(i,j,k) != a->subgrid(i,j,k+1);
        isupper = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j,k+1);
        islower  = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j,k+1);
        isleft  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i,j,k+1);
        isright = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i,j,k+1);

        condperiod = (p->C23 == 2) && (a->flag(i,j,k+1) == -21);
        isupperperiod = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);
        isleftperiod  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);
        isrightperiod = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);

        if((condition && islower && isleft)||(condperiod && islowerperiod && isleftperiod))
        paracell[0][0][2]+=1;

        if((condition && isupper && isleft)||(condperiod && isupperperiod && isleftperiod))
        paracell[0][2][2]+=1;

        if((condition && islower && isright)||(condperiod && islowerperiod && isrightperiod))
        paracell[2][0][2]+=1;

        if((condition && isupper && isright)||(condperiod && isupperperiod && isrightperiod))
        paracell[2][2][2]+=1;


        for(n=0;n<3;++n)
        for(q=0;q<3;++q)
        for(r=0;r<3;++r)
        if(paracell[n][q][r]==3)
        paracell[n][q][r]=3;

        //---------------------------------------------------------------------------------------------
        //  1
        condition =  a->flag(i,j,k)>0 && a->flag(i-1,j,k)>0 && a->subgrid(i,j,k) != a->subgrid(i-1,j,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i-1,j,k);
        isright   = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i-1,j,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i-1,j,k);
        isleft    = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i-1,j,k);

        condperiod = (p->C21 == 2) && (a->flag(i-1,j,k) == -21);
        isupperperiod   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isrightperiod    = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod    = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleftperiod     = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);

        if((condition && isright && islower)|(condperiod && isrightperiod && islowerperiod))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j+1;
            a->para1co[a->paraco1count][2]=k-1;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=2;
            a->para1co[a->paraco1count][5]=5;
            a->para1co[a->paraco1count][6]=paracell[0][2][0];
            a->paraco1count++;
        }

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j+1;
            a->para1co[a->paraco1count][2]=k+1;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=2;
            a->para1co[a->paraco1count][5]=6;
            a->para1co[a->paraco1count][6]=paracell[0][2][2];
            a->paraco1count++;
        }

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j-1;
            a->para1co[a->paraco1count][2]=k-1;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=3;
            a->para1co[a->paraco1count][5]=5;
            a->para1co[a->paraco1count][6]=paracell[0][0][0];
            a->paraco1count++;
        }

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        {
            a->para1co[a->paraco1count][0]=i;
            a->para1co[a->paraco1count][1]=j-1;
            a->para1co[a->paraco1count][2]=k+1;
            a->para1co[a->paraco1count][3]=a->subgrid(i,j,k);
            a->para1co[a->paraco1count][4]=3;
            a->para1co[a->paraco1count][5]=6;
            a->para1co[a->paraco1count][6]=paracell[0][0][2];
            a->paraco1count++;
        }

        //  4
        condition = a->flag(i,j,k)>0 && a->flag(i+1,j,k)>0 && a->subgrid(i,j,k) != a->subgrid(i+1,j,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i+1,j,k);
        isright   = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i+1,j,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i+1,j,k);
        isleft    = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i+1,j,k);

        condperiod = (p->C21 == 2) && (a->flag(i+1,j,k) == -21);
        isupperperiod   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isrightperiod   = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleftperiod    = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j+1;
            a->para4co[a->paraco4count][2]=k-1;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=2;
            a->para4co[a->paraco4count][5]=5;
            a->para4co[a->paraco4count][6]=paracell[2][2][0];
            a->paraco4count++;
        }

        if((condition && isright && isupper)||(condperiod && isrightperiod && islowerperiod))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j+1;
            a->para4co[a->paraco4count][2]=k+1;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=2;
            a->para4co[a->paraco4count][5]=6;
            a->para4co[a->paraco4count][6]=paracell[2][2][2];
            a->paraco4count++;
        }

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j-1;
            a->para4co[a->paraco4count][2]=k-1;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=3;
            a->para4co[a->paraco4count][5]=5;
            a->para4co[a->paraco4count][6]=paracell[2][0][0];
            a->paraco4count++;
        }

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        {
            a->para4co[a->paraco4count][0]=i;
            a->para4co[a->paraco4count][1]=j-1;
            a->para4co[a->paraco4count][2]=k+1;
            a->para4co[a->paraco4count][3]=a->subgrid(i,j,k);
            a->para4co[a->paraco4count][4]=3;
            a->para4co[a->paraco4count][5]=6;
            a->para4co[a->paraco4count][6]=paracell[2][0][2];
            a->paraco4count++;
        }

        //  3
        condition = a->flag(i,j,k)>0 && a->flag(i,j-1,k)>0 && a->subgrid(i,j,k) != a->subgrid(i,j-1,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j-1,k);
        isright   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j-1,k)!=a->subgrid(i,j-1,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j-1,k);
        isleft    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j-1,k)!=a->subgrid(i,j-1,k);

        condperiod = (p->C22 == 2) && (a->flag(i,j-1,k) == -21);
        isupperperiod   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isrightperiod   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);
        islowerperiod   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleftperiod    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        {
            a->para3co[a->paraco3count][0]=i-1;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k-1;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=1;
            a->para3co[a->paraco3count][5]=5;
            a->para3co[a->paraco3count][6]=paracell[0][0][0];
            a->paraco3count++;
        }

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        {
            a->para3co[a->paraco3count][0]=i-1;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k+1;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=1;
            a->para3co[a->paraco3count][5]=6;
            a->para3co[a->paraco3count][6]=paracell[0][0][2];
            a->paraco3count++;
        }

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        {
            a->para3co[a->paraco3count][0]=i+1;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k-1;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=4;
            a->para3co[a->paraco3count][5]=5;
            a->para3co[a->paraco3count][6]=paracell[2][0][0];
            a->paraco3count++;
        }

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        {
            a->para3co[a->paraco3count][0]=i+1;
            a->para3co[a->paraco3count][1]=j;
            a->para3co[a->paraco3count][2]=k+1;
            a->para3co[a->paraco3count][3]=a->subgrid(i,j,k);
            a->para3co[a->paraco3count][4]=4;
            a->para3co[a->paraco3count][5]=6;
            a->para3co[a->paraco3count][6]=paracell[2][0][2];
            a->paraco3count++;
        }

        //  2
        condition = a->flag(i,j,k)>0 && a->flag(i,j+1,k)>0 && a->subgrid(i,j,k) != a->subgrid(i,j+1,k);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j+1,k);
        isright   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j+1,k)!=a->subgrid(i,j+1,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j+1,k);
        isleft    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j+1,k)!=a->subgrid(i,j+1,k);

        condperiod = (p->C22 == 2) && (a->flag(i,j+1,k) == -21);
        isupper   = a->subgrid(i,j,k+1)!=a->subgrid(i,j,k);
        isright   = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);
        islower   = a->subgrid(i,j,k-1)!=a->subgrid(i,j,k);
        isleft    = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);

        if((condition && isleft && islower)||(condperiod && isleftperiod && islowerperiod))
        {
            a->para2co[a->paraco2count][0]=i-1;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k-1;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=1;
            a->para2co[a->paraco2count][5]=5;
            a->para2co[a->paraco2count][6]=paracell[0][2][0];
            a->paraco2count++;
        }

        if((condition && isleft && isupper)||(condperiod && isleftperiod && isupperperiod))
        {
            a->para2co[a->paraco2count][0]=i-1;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k+1;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=1;
            a->para2co[a->paraco2count][5]=6;
            a->para2co[a->paraco2count][6]=paracell[0][2][2];
            a->paraco2count++;
        }

        if((condition && isright && islower)||(condperiod && isrightperiod && islowerperiod))
        {
            a->para2co[a->paraco2count][0]=i+1;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k-1;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=4;
            a->para2co[a->paraco2count][5]=5;
            a->para2co[a->paraco2count][6]=paracell[2][2][0];
            a->paraco2count++;
        }

        if((condition && isright && isupper)||(condperiod && isrightperiod && isupperperiod))
        {
            a->para2co[a->paraco2count][0]=i+1;
            a->para2co[a->paraco2count][1]=j;
            a->para2co[a->paraco2count][2]=k+1;
            a->para2co[a->paraco2count][3]=a->subgrid(i,j,k);
            a->para2co[a->paraco2count][4]=4;
            a->para2co[a->paraco2count][5]=6;
            a->para2co[a->paraco2count][6]=paracell[2][2][2];
            a->paraco2count++;
        }

        //  5
        condition = a->flag(i,j,k)>0 && a->flag(i,j,k-1)>0 && a->subgrid(i,j,k) != a->subgrid(i,j,k-1);
        isupper = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k-1)!=a->subgrid(i,j,k-1);
        islower  = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k-1)!=a->subgrid(i,j,k-1);
        isleft  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k-1)!=a->subgrid(i,j,k-1);
        isright = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k-1)!=a->subgrid(i,j,k-1);

        condperiod = (p->C23 == 2) && (a->flag(i,j,k-1) == -21);
        isupperperiod = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);
        isleftperiod  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);
        isrightperiod = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);

        if((condition && islower && isleft)||(condperiod && islowerperiod && isleftperiod))
        {
            a->para5co[a->paraco5count][0]=i-1;
            a->para5co[a->paraco5count][1]=j-1;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=1;
            a->para5co[a->paraco5count][5]=3;
            a->para5co[a->paraco5count][6]=paracell[0][0][0];
            a->paraco5count++;
        }

        if((condition && isupper && isleft)||(condperiod && isupperperiod && isleftperiod))
        {
            a->para5co[a->paraco5count][0]=i-1;
            a->para5co[a->paraco5count][1]=j+1;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=1;
            a->para5co[a->paraco5count][5]=2;
            a->para5co[a->paraco5count][6]=paracell[0][2][0];
            a->paraco5count++;
        }

        if((condition && islower && isright)||(condperiod && islowerperiod && isrightperiod))
        {
            a->para5co[a->paraco5count][0]=i+1;
            a->para5co[a->paraco5count][1]=j-1;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=4;
            a->para5co[a->paraco5count][5]=3;
            a->para5co[a->paraco5count][6]=paracell[2][0][0];
            a->paraco5count++;
        }

        if((condition && isupper && isright)||(condperiod && isupperperiod && isrightperiod))
        {
            a->para5co[a->paraco5count][0]=i+1;
            a->para5co[a->paraco5count][1]=j+1;
            a->para5co[a->paraco5count][2]=k;
            a->para5co[a->paraco5count][3]=a->subgrid(i,j,k);
            a->para5co[a->paraco5count][4]=4;
            a->para5co[a->paraco5count][5]=2;
            a->para5co[a->paraco5count][6]=paracell[2][2][0];
            a->paraco5count++;
        }

        //  6
        condition = a->flag(i,j,k)>0 && a->flag(i,j,k+1)>0 && a->subgrid(i,j,k) != a->subgrid(i,j,k+1);
        isupper = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j+1,k+1)!=a->subgrid(i,j,k+1);
        islower  = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k)  && a->subgrid(i,j-1,k+1)!=a->subgrid(i,j,k+1);
        isleft  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i-1,j,k+1)!=a->subgrid(i,j,k+1);
        isright = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k)  && a->subgrid(i+1,j,k+1)!=a->subgrid(i,j,k+1);

        condperiod = (p->C23 == 2) && (a->flag(i,j,k+1) == -21);
        isupperperiod = a->subgrid(i,j+1,k)!=a->subgrid(i,j,k);
        islowerperiod = a->subgrid(i,j-1,k)!=a->subgrid(i,j,k);
        isleftperiod  = a->subgrid(i-1,j,k)!=a->subgrid(i,j,k);
        isrightperiod = a->subgrid(i+1,j,k)!=a->subgrid(i,j,k);

        if((condition && islower && isleft)||(condperiod && islowerperiod && isleftperiod))
        {
            a->para6co[a->paraco6count][0]=i-1;
            a->para6co[a->paraco6count][1]=j-1;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=1;
            a->para6co[a->paraco6count][5]=3;
            a->para6co[a->paraco6count][6]=paracell[0][0][2];
            a->paraco6count++;
        }

        if((condition && isupper && isleft)||(condperiod && isupperperiod && isleftperiod))
        {
            a->para6co[a->paraco6count][0]=i-1;
            a->para6co[a->paraco6count][1]=j+1;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=1;
            a->para6co[a->paraco6count][5]=2;
            a->para6co[a->paraco6count][6]=paracell[0][2][2];
            a->paraco6count++;
        }

        if((condition && islower && isright)||(condperiod && islowerperiod && isrightperiod))
        {
            a->para6co[a->paraco6count][0]=i+1;
            a->para6co[a->paraco6count][1]=j-1;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=4;
            a->para6co[a->paraco6count][5]=3;
            a->para6co[a->paraco6count][6]=paracell[2][0][2];
            a->paraco6count++;
        }

        if((condition && isupper && isright)||(condperiod && isupperperiod && isrightperiod))
        {
            a->para6co[a->paraco6count][0]=i+1;
            a->para6co[a->paraco6count][1]=j+1;
            a->para6co[a->paraco6count][2]=k;
            a->para6co[a->paraco6count][3]=a->subgrid(i,j,k);
            a->para6co[a->paraco6count][4]=4;
            a->para6co[a->paraco6count][5]=2;
            a->para6co[a->paraco6count][6]=paracell[2][2][2];
            a->paraco6count++;
        }
    } // END LOOP

    ddout<<"co1: "<<a->paraco1count<<" co2: "<<a->paraco2count<<" co3: "<<a->paraco3count<<" co4: "<<a->paraco4count<<" co5: "<<a->paraco5count<<" co6: "<<a->paraco6count<<endl;
}
}
