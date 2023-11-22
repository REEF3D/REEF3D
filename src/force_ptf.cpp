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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"force_ptf.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"iterators1D.h"

force_ptf::force_ptf(lexer* p, fdm_ptf *e, ghostcell *pgc,int(qn)):ID(qn){}

force_ptf::~force_ptf(){}

void force_ptf::ini(lexer *p, fdm_ptf *e, ghostcell *pgc)
{
    force_ptfprintcount=0;
    
    print_ini_ptf(p,e,pgc);
}

void force_ptf::start(lexer *p, fdm_ptf *e, ghostcell *pgc)
{
    F_x_tot=0.0;
    F_y_tot=0.0;
    F_z_tot=0.0;
    
   OSOLIDLOOP
   {
        if(p->flag4[Im1JK]>0)
           F_x_tot+=e->press(i-1,j,k)*p->DYN[JP]*p->DZN[KP];
    
        if(p->flag4[Ip1JK]>0)
           F_x_tot-=e->press(i+1,j,k)*p->DYN[JP]*p->DZN[KP];
           
        if(p->flag4[IJm1K]>0)
           F_y_tot+=e->press(i,j-1,k)*p->DXN[IP]*p->DZN[KP];
           
        if(p->flag4[IJp1K]>0)
           F_y_tot-=e->press(i,j+1,k)*p->DXN[IP]*p->DZN[KP];
           
        if(p->flag4[IJKm1]>0)
            F_z_tot+=e->press(i,j,k-1)*p->DXN[IP]*p->DYN[JP];
            
        if(p->flag4[IJKp1]>0)
            F_z_tot-=e->press(i,j,k+1)*p->DXN[IP]*p->DYN[JP];
   }
    // Sum up to distribute forces
    F_x_tot = pgc->globalsum(F_x_tot);
    F_y_tot = pgc->globalsum(F_y_tot);
    F_z_tot = pgc->globalsum(F_z_tot);

    // Print
    if(p->mpirank==0)
    {
        print_force_ptf(p,e,pgc);
    }
}


