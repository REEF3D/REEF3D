/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"net_interface.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>

#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"

void net_interface::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if (p->X320 == 0)
    pnet.push_back(new net_void());
    
    else
    {
        pgc->bcast_int(&p->net_count, 1);
        
        if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_6DOF_Net",0777);	

        else
        p->X320_type = new int[p->net_count];

		
        pnet.reserve(p->net_count);	
  
		NETLOOP
		{
            pgc->bcast_int(&p->X320_type[n], 1);
			
            if(p->X320_type[n]>10)
            pnet.push_back(new net_barDyn(n,p));

            else if(p->X320_type[n]<4)
            pnet.push_back(new net_barQuasiStatic(n,p));

            else
            pnet.push_back(new net_sheet(n,p));

            pnet[n]->initialize_nhflow(p,d,pgc);
		}
    }

}