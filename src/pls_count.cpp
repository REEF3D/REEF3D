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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_pls::parcount(lexer* p, fdm* a, ghostcell* pgc)
{		
	LOOP
	{
		posnum(i,j,k)=0.0;
		negnum(i,j,k)=0.0;
	}
	
	pgc->start4(p,posnum,1);
	pgc->start4(p,negnum,1);
		
	// POS
    for(n=0;n<posactive;++n)
    if(posflag[n]>0)
    {
        i=int((pos[n][0])/dx);
        j=int((pos[n][1])/dx);
        k=int((pos[n][2])/dx);

    posnum(i,j,k)+=1.0;
    }

    // NEG
    for(n=0;n<negactive;++n)
    if(negflag[n]>0)
    {
        i=int((neg[n][0])/dx);
        j=int((neg[n][1])/dx);
        k=int((neg[n][2])/dx);

    negnum(i,j,k)+=1.0;
    }

}
