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
#include<math.h>

void particle_pls::remove(lexer* p, fdm* a, ghostcell* pgc)
{
    removed = 0;

    for(n=0;n<posactive;++n)
    {
        // POS
        if(posflag[n]>0)
        {
            i=int((pos[n][0])/dx);
            j=int((pos[n][1])/dx);
            k=int((pos[n][2])/dx);

            check=boundcheck(p,a,i,j,k,1);
			
			// remove particle_plss too far away from ls
            if(check==1)
            if(p->flag5[IJK]>0  || fabs(pos[n][3])>epsi)
            {
			pcount++;
            posflag[n]=0;
            posmem[pcount]=n;
            removed++;
            }
			
			// remove out of bounds particle_plss
            if(check==0)
            {
			pcount++;
            posflag[n]=0;
            posmem[pcount]=n;
            removed++;
            }
        }
    }

    for(n=0;n<negactive;++n)
    {

        //NEG
        if(negflag[n]>0)
        {
            i=int((neg[n][0])/dx);
            j=int((neg[n][1])/dx);
            k=int((neg[n][2])/dx);

            check=boundcheck(p,a,i,j,k,1);

            if(check==1)
            if(p->flag5[IJK]>0 || fabs(neg[n][3])>epsi)
            {
			ncount++;
            negflag[n]=0;
            negmem[ncount]=n;
            removed++;
            }

            if(check==0)
            {
			ncount++;
            negflag[n]=0;
            negmem[ncount]=n;
            removed++;
            }
        }
    }
}
