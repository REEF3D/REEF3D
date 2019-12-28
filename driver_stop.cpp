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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"printer.h"

void driver::stop(lexer *p, fdm *a, ghostcell *pgc)
{	
    if(p->umax>p->N61 || p->vmax>p->N61 || p->wmax>p->N61)
    {
    
        if(p->mpirank==0)
        cout<<endl<<"EMERGENCY STOP  --  velocities exceeding critical value N 61"<<endl<<endl;
        
    
     if(p->A10==5)
     pprint->print_vtu(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,psed);
     
     pgc->final();
     exit(0);
    }
    
}