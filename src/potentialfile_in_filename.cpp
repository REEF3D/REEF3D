
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

#include"potentialfile_in.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void potentialfile_in::filename(lexer *p, fdm *a, ghostcell *pgc, char *name, int num)
{
	
		if(p->gcin_count>0)
		{
			// open file
			if(p->P14==0)
			{
			if(num<10)
			sprintf(name,"REEF3D-flowfile-%i-00000%i.r3d",p->I230,num);
 
			if(num<100&&num>9)
			sprintf(name,"REEF3D-flowfile-%i-0000%i.r3d",p->I230,num);

			if(num<1000&&num>99)
			sprintf(name,"REEF3D-flowfile-%i-000%i.r3d",p->I230,num);

			if(num<10000&&num>999)
			sprintf(name,"REEF3D-flowfile-%i-00%i.r3d",p->I230,num);

			if(num<100000&&num>9999)
			sprintf(name,"REEF3D-flowfile-%i-0%i.r3d",p->I230,num);

			if(num>99999)
			sprintf(name,"REEF3D-flowfile-%i-%i.r3d",p->I230,num);
			}
			
			if(p->P14==1)
			{
			if(num<10)
			sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%i-00000%i.r3d",p->I230,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%i-0000%i.r3d",p->I230,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%i-000%i.r3d",p->I230,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%i-00%i.r3d",p->I230,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%i-0%i.r3d",p->I230,num);

			if(num>99999)
			sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%i-%i.r3d",p->I230,num);
			}
		}
    
}


