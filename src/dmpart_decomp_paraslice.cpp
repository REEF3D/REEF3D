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

void decomp::paraslicesurface(lexer* p, dive* a)
{
    a->paraslice1count=0;
    a->paraslice2count=0;
    a->paraslice3count=0;
    a->paraslice4count=0;

    XYLOOP
    {
        if(a->subslice(i-1,j)!=a->subslice(i,j) && a->flagslice(i-1,j)>0 && a->flagslice(i,j)>0)
        {
            a->paraslice1sf[a->paraslice1count][0]=i;
            a->paraslice1sf[a->paraslice1count][1]=j;
            a->paraslice1count++;
        }

        if(a->subslice(i+1,j)!=a->subslice(i,j) && a->flagslice(i+1,j)>0 && a->flagslice(i,j)>0)
        {
            a->paraslice4sf[a->paraslice4count][0]=i;
            a->paraslice4sf[a->paraslice4count][1]=j;
            a->paraslice4count++;
        }

        if(a->subslice(i,j-1)!=a->subslice(i,j) && a->flagslice(i,j-1)>0 && a->flagslice(i,j)>0)
        {
            a->paraslice3sf[a->paraslice3count][0]=i;
            a->paraslice3sf[a->paraslice3count][1]=j;
            a->paraslice3count++;
        }

        if(a->subslice(i,j+1)!=a->subslice(i,j) && a->flagslice(i,j+1)>0 && a->flagslice(i,j)>0)
        {
            a->paraslice2sf[a->paraslice2count][0]=i;
            a->paraslice2sf[a->paraslice2count][1]=j;
            a->paraslice2count++;
        }
    }
}

void decomp::paraslicecosurface(lexer* p, dive* a)
{
    a->paracoslice1count=0;
    a->paracoslice2count=0;
    a->paracoslice3count=0;
    a->paracoslice4count=0;

    XYLOOP
    {
        //-----------------------------------------------------------------------------------------------
        //  1
        if(a->subslice(i-1,j)!=a->subslice(i,j)
        && a->flagslice(i,j)>0 && a->flagslice(i-1,j)>0
        && ((a->flagslice(i,j-1)<0  && a->flagslice(i-1,j-1)<0)
            ||(a->subslice(i,j-1)!=a->subslice(i,j)  && a->subslice(i-1,j-1)!=a->subslice(i-1,j)) ))
        {
        a->paracoslice1sf[a->paracoslice1count][0]=i;
        a->paracoslice1sf[a->paracoslice1count][1]=j-1;
        a->paracoslice1sf[a->paracoslice1count][2]=a->subslice(i,j);
        a->paracoslice1sf[a->paracoslice1count][3]=3;
        ++a->paracoslice1count;
        }

        if(a->subslice(i-1,j)!=a->subslice(i,j)
        && a->flagslice(i,j)>0
        && a->flagslice(i-1,j)>0
        && ((a->flagslice(i,j+1)<0  && a->flagslice(i-1,j+1)<0)
            ||(a->subslice(i,j+1)!=a->subslice(i,j)  && a->subslice(i-1,j+1)!=a->subslice(i-1,j))))
        {
        a->paracoslice1sf[a->paracoslice1count][0]=i;
        a->paracoslice1sf[a->paracoslice1count][1]=j+1;
        a->paracoslice1sf[a->paracoslice1count][2]=a->subslice(i,j);
        a->paracoslice1sf[a->paracoslice1count][3]=2;
        a->paracoslice1count++;
        }

        //  4
        if(a->subslice(i+1,j)!=a->subslice(i,j)
        && a->flagslice(i,j)>0
        && a->flagslice(i+1,j)>0
        && ((a->flagslice(i,j-1)<0  && a->flagslice(i+1,j-1)<0)
            ||(a->subslice(i,j-1)!=a->subslice(i,j)  && a->subslice(i+1,j-1)!=a->subslice(i+1,j))))
        {
            a->paracoslice4sf[a->paracoslice4count][0]=i;
            a->paracoslice4sf[a->paracoslice4count][1]=j-1;
            a->paracoslice4sf[a->paracoslice4count][2]=a->subslice(i,j);
            a->paracoslice4sf[a->paracoslice4count][3]=3;
            a->paracoslice4count++;
        }


        if(a->subslice(i+1,j)!=a->subslice(i,j)
        && a->flagslice(i,j)>0
        && a->flagslice(i+1,j)>0
        && ((a->flagslice(i,j+1)<0  && a->flagslice(i+1,j+1)<0)
            ||(a->subslice(i,j+1)!=a->subslice(i,j)  && a->subslice(i+1,j+1)!=a->subslice(i+1,j))))
        {
            a->paracoslice4sf[a->paracoslice4count][0]=i;
            a->paracoslice4sf[a->paracoslice4count][1]=j+1;
            a->paracoslice4sf[a->paracoslice4count][2]=a->subslice(i,j);
            a->paracoslice4sf[a->paracoslice4count][3]=2;
            a->paracoslice4count++;
        }

        //  3
        if(a->subslice(i,j-1)!=a->subslice(i,j)
        && a->flagslice(i,j)>0
        && a->flagslice(i,j-1)>0
        && ((a->flagslice(i-1,j)<0  && a->flagslice(i-1,j-1)<0)
            ||(a->subslice(i-1,j)!=a->subslice(i,j)  && a->subslice(i-1,j-1)!=a->subslice(i,j-1))))
        {
            a->paracoslice3sf[a->paracoslice3count][0]=i-1;
            a->paracoslice3sf[a->paracoslice3count][1]=j;
            a->paracoslice3sf[a->paracoslice3count][2]=a->subslice(i,j);
            a->paracoslice3sf[a->paracoslice3count][3]=1;
            a->paracoslice3count++;
        }

        if(a->subslice(i,j-1)!=a->subslice(i,j)
        && a->flagslice(i,j)>0
        && a->flagslice(i,j-1)>0
        && ((a->flagslice(i+1,j)<0  && a->flagslice(i+1,j-1)<0)
            ||(a->subslice(i+1,j)!=a->subslice(i,j)  && a->subslice(i+1,j-1)!=a->subslice(i,j-1))))
        {
            a->paracoslice3sf[a->paracoslice3count][0]=i+1;
            a->paracoslice3sf[a->paracoslice3count][1]=j;
            a->paracoslice3sf[a->paracoslice3count][2]=a->subslice(i,j);
            a->paracoslice3sf[a->paracoslice3count][3]=4;
            a->paracoslice3count++;
        }

        //  2
        if(a->subslice(i,j+1)!=a->subslice(i,j)

        && a->flagslice(i,j)>0
        && a->flagslice(i,j+1)>0
        && ((a->flagslice(i-1,j)<0  && a->flagslice(i-1,j+1)<0)
            ||(a->subslice(i-1,j)!=a->subslice(i,j)  && a->subslice(i-1,j+1)!=a->subslice(i,j+1))))
        {
            a->paracoslice2sf[a->paracoslice2count][0]=i-1;
            a->paracoslice2sf[a->paracoslice2count][1]=j;
            a->paracoslice2sf[a->paracoslice2count][2]=a->subslice(i,j);
            a->paracoslice2sf[a->paracoslice2count][3]=1;
            a->paracoslice2count++;
        }

        if(a->subslice(i,j+1)!=a->subslice(i,j)
        && a->flagslice(i,j)>0
        && a->flagslice(i,j+1)>0
        && ((a->flagslice(i+1,j)<0  && a->flagslice(i+1,j+1)<0)
            ||(a->subslice(i+1,j)!=a->subslice(i,j)  && a->subslice(i+1,j+1)!=a->subslice(i,j+1))))
        {
            a->paracoslice2sf[a->paracoslice2count][0]=i+1;
            a->paracoslice2sf[a->paracoslice2count][1]=j;
            a->paracoslice2sf[a->paracoslice2count][2]=a->subslice(i,j);
            a->paracoslice2sf[a->paracoslice2count][3]=4;
            a->paracoslice2count++;
        }
    }
}
}
