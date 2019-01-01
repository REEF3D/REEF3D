/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include"slice.h"

void ghostcell::dgcslpol(lexer* p,slice& f,int **dgc, int dgc_count, int gcv)
{
    double val0,val1,val2,val3,val4;
    int aa,bb,cc;
    int acheck,bcheck;
    int a1,a2,b1,b2;

    val0=val1=val2=val3=val4=0.0;

    //cout<<p->mpirank<<" DGC: "<<dgc_count<<" gcv: "<<gcv<<endl;

    for(n=0;n<dgc_count;++n)
    {
    i=dgc[n][0];
    j=dgc[n][1];


    acheck=bcheck=0;

    // pre-check
    for(q=0;q<dgc[n][3];++q)
    {
        if(dgc[n][4+q]==1)
        ++acheck;

        if(dgc[n][4+q]==4)
        ++acheck;

        if(dgc[n][4+q]==2)
        ++bcheck;

        if(dgc[n][4+q]==3)
        ++bcheck;
    }

// -------------------------------
    if(dgc[n][3]==2 && acheck!=2 && bcheck!=2)
    {
        aa=bb=cc=0;
        val1=val2=val3=val4=0.0;

        for(q=0;q<2;++q)
        {
            if(dgc[n][4+q]==1)
            {
            val1=f(i-1,j);
            aa=-1;
            }

            if(dgc[n][4+q]==2)
            {
            val2=f(i,j+1);
            bb=1;
            }

            if(dgc[n][4+q]==3)
            {
            val3=f(i,j-1);
            bb=-1;
            }

            if(dgc[n][4+q]==4)
            {
            val4=f(i+1,j);
            aa=1;
            }
        }
        f(i+aa,j+bb)=0.5*(val1+val2+val3+val4);

    }
    }
}
