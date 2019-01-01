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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::surfarea(lexer* p, fdm *a, ghostcell *pgc)
{
    double A=0.0;
	int fnumg=0;

    for(n=0;n<fnum;++n)
    {
        if(surfnum[n]==3)
        farea[n] = area(0,1,2);

        if(surfnum[n]==4)
        farea[n] = area(0,1,2) + area(0,2,3);

        if(surfnum[n]==5)
        farea[n] = area(0,1,4) + area(1,2,3) + area(1,3,4);

        if(surfnum[n]==6)
        farea[n] = area(0,1,5) + area(1,2,5) + area(2,3,4) + area(2,4,5);

        A+=farea[n];
    }
    A  = pgc->globalsum(A);
	fnumg = pgc->globalisum(fnum);

    if(p->mpirank==0)
    cout<<"Atot pier:  "<<A<<"  fnum_obj: "<<fnumg<<endl;
}

double force::area(int q1, int q2, int q3)
{
    double area=0.0;
    double s,a,b,c;

    a = sqrt(pow(ccpt[n][q2][0]-ccpt[n][q1][0], 2.0) + pow(ccpt[n][q2][1]-ccpt[n][q1][1], 2.0) + pow(ccpt[n][q2][2]-ccpt[n][q1][2], 2.0));

    b = sqrt(pow(ccpt[n][q3][0]-ccpt[n][q2][0], 2.0) + pow(ccpt[n][q3][1]-ccpt[n][q2][1], 2.0) + pow(ccpt[n][q3][2]-ccpt[n][q2][2], 2.0));

    c = sqrt(pow(ccpt[n][q1][0]-ccpt[n][q3][0], 2.0) + pow(ccpt[n][q1][1]-ccpt[n][q3][1], 2.0) + pow(ccpt[n][q1][2]-ccpt[n][q3][2], 2.0));

    s = 0.5*(a+b+c);

    area = sqrt(s*(s-a)*(s-b)*(s-c));

    return area;
}




