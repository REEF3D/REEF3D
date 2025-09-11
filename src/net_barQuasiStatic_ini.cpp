/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2025 Tobias Martin

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"net_barQuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"	

void net_barQuasiStatic::initialize_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    //- Initialise net model
    if (p->X320_type[nNet]==1)
    {
        bag_ini(p,pgc);
        
        buildNet_bag(p);
    }
    else if (p->X320_type[nNet]==2)   
    {
        cyl_ini(p,pgc);
        
        buildNet_cyl(p); 
    }
    else if (p->X320_type[nNet]==3)   
    {
        wall_ini(p,pgc);
        
        buildNet_wall(p);    
    }  
    
    //- Update porous zone
    coupling_dlm_cfd(p,a,pgc);

    print(p);
}

void net_barQuasiStatic::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    //- Initialise net model
    if (p->X320_type[nNet]==1)
    {
        bag_ini(p,pgc);
        
        buildNet_bag(p);
    }
    else if (p->X320_type[nNet]==2)   
    {
        cyl_ini(p,pgc);
        
        buildNet_cyl(p); 
    }
    else if (p->X320_type[nNet]==3)   
    {
        wall_ini(p,pgc);
        
        buildNet_wall(p);    
    }  
    
    //- Update porous zone
    coupling_dlm_nhflow(p,d,pgc);

    print(p);
}
