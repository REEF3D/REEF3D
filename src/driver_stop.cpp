/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include"nhflow_printer.h"
#include"fnpf_printer.h"

void driver::stop(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->A10==4 || p->A10==6)
    {
        bool check=false;
        
        ULOOP
            if(a->u(i,j,k)!=a->u(i,j,k))
                check=true;
        
        VLOOP
            if(a->v(i,j,k)!=a->v(i,j,k))
                check=true;
        
        WLOOP
            if(a->w(i,j,k)!=a->w(i,j,k))
                check=true;
        
        LOOP
            if(a->press(i,j,k)!=a->press(i,j,k))
                check=true;
        
        if(check)
        {
            if(p->mpirank==0)
                cout<<endl<<"EMERGENCY STOP  --  solver breaking down - NAN values\n"<<endl;

            pprint->print_stop(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);
            
            pgc->final();
            exit(1);
        }
    }
        
    if(p->umax>p->N61 || p->vmax>p->N61 || p->wmax>p->N61)
    {
        if(p->mpirank==0)
            cout<<endl<<"EMERGENCY STOP  --  velocities exceeding critical value N 61\n"<<endl;
    
        if(p->A10==3)
            pfprint->print_stop(p,c,pgc);
    
        if(p->A10==4 || p->A10==6)
            pprint->print_stop(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);
        
        if(p->A10==5)
            pnhfprint->print_stop(p,d,pgc,pflow,pnhfturb,psed);
     
        pgc->final();
        exit(1);
    }
    
    // Solver Status
    p->solver_status = pgc->globalimax(p->solver_status);
    /*
    if(p->solver_status>=1)
    {
        if(p->mpirank==0)
            cout<<endl<<" HYPRE solver broke down! Emergency Stop! "<<p->solver_status<<"\n"<<endl;
        
        if(p->A10==3)
            pfprint->print_vtu(p,c,pgc);
        
        if(p->A10==4 || p->A10==5 || p->A10==6)
            pprint->print_vtu(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);
        
        pgc->final();
        exit(1);
    }
    */
}
