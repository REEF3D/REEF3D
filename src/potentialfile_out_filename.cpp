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

#include"potentialfile_out.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void potentialfile_out::filename(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
	int num;

    num = p->count;

			// open file
			if(p->P14==0)
			{
			if(num<10)
			sprintf(name,"REEF3D-potentialfile-%i-00000%i.r3d",n+1,num);
 
			if(num<100&&num>9)
			sprintf(name,"REEF3D-potentialfile-%i-0000%i.r3d",n+1,num);

			if(num<1000&&num>99)
			sprintf(name,"REEF3D-potentialfile-%i-000%i.r3d",n+1,num);

			if(num<10000&&num>999)
			sprintf(name,"REEF3D-potentialfile-%i-00%i.r3d",n+1,num);

			if(num<100000&&num>9999)
			sprintf(name,"REEF3D-potentialfile-%i-0%i.r3d",n+1,num);

			if(num>99999)
			sprintf(name,"REEF3D-potentialfile-%i-%i.r3d",n+1,num);
			}
			
			if(p->P14==1)
			{
			if(num<10)
			sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialfile-%i-00000%i.r3d",n+1,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialfile-%i-0000%i.r3d",n+1,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialfile-%i-000%i.r3d",n+1,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialfile-%i-00%i.r3d",n+1,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialfile-%i-0%i.r3d",n+1,num);

			if(num>99999)
			sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialfile-%i-%i.r3d",n+1,num);
			}
		
    
}


