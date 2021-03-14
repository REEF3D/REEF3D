/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

void force::surfarea_gcb(lexer* p, fdm *a, ghostcell *pgc)
{
    double A=0.0;
	int fnumg=0;

    for(n=0;n<fnum;++n)
    {
        A += farea[n];
    }
    
    A  = pgc->globalsum(A);
	fnumg = pgc->globalisum(fnum);

    if(p->mpirank==0)
    cout<<"Atot pier:  "<<A<<"  fnum: "<<fnumg<<endl;
}





