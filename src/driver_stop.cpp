/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"printer.h"
#include"waves_header.h"
#include"nhflow_header.h"

void driver::stop(lexer *p, fdm *a, ghostcell *pgc)
{	 
    
    if(p->A10==4 || p->A10==6)
    {
    int check=0;
    
    ULOOP
    if(a->u(i,j,k)!=a->u(i,j,k))
    check=1;
    
    VLOOP
    if(a->v(i,j,k)!=a->v(i,j,k))
    check=1;
    
    WLOOP
    if(a->w(i,j,k)!=a->w(i,j,k))
    check=1;
    
    LOOP
    if(a->press(i,j,k)!=a->press(i,j,k))
    check=1;
    
    if(check==1)
    {
    
        if(p->mpirank==0)
        cout<<endl<<"!!!  EMERGENCY STOP  --  solver breaking down - NAN values  !!!"<<endl<<endl;

     pprint->print_stop(p,a,pgc,pturb,pheat,pflow,pdata,pconc,pmp,psed);
     
     //pgc->final(true);
    }
    }
        
    if(p->umax>p->N61 || p->vmax>p->N61 || p->wmax>p->N61)
    {
    
        if(p->mpirank==0)
        {
        cout<<endl<<"!!! EMERGENCY STOP  --  velocities exceeding critical value N 61  !!!"<<endl<<endl;
        
        cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
        cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
        cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
        }
    
        if(p->A10==3)
        pprint->print_stop(p,c,pgc);
        
        if(p->A10==4 || p->A10==6)
        pprint->print_stop(p,a,pgc,pturb,pheat,pflow,pdata,pconc,pmp,psed);
        
        if(p->A10==5)
        pprint->print_stop(p,d,pgc,pflow,pnhfturb,psed);
        
        pgc->final(true);
    }
    
    // Solver Status
    p->solver_status = pgc->globalimax(p->solver_status);
    p->solver_error = pgc->globalimax(p->solver_error);
    
    if(p->solver_error>=1)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! EMERGENCY STOP  --  HYPRE solver broke down!  !!!     "<<p->solver_error<<endl<<endl;
        
        if(p->A10==3)
        pprint->print_stop(p,c,pgc);
        
        if(p->A10==4 || p->A10==5 || p->A10==6)
        pprint->print_stop(p,a,pgc,pturb,pheat,pflow,pdata,pconc,pmp,psed);

        if(p->A10==5)
        pprint->print_stop(p,d,pgc,pflow,pnhfturb,psed);
        
        pgc->final(true);
    }
}
