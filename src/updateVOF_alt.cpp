/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::updateVOF_alt(fdm* a, lexer* p)
{
    double V_w_celltotal, V_a_celltotal;
	LOOP
    {
        V_w_celltotal=V_w_old(i,j,k)+V_w_update(i,j,k);
        V_a_celltotal=V_a_old(i,j,k)+V_a_update(i,j,k);
        if(V_w_celltotal<0.0)
        {
            cout<<"water vol loss"<<endl;
            a->vof(i,j,k)=0.0;
        }
        /*else if(Vc_a_celltotal<0.0)
        {
            cout<<"air vol loss"<<endl;
            a->vof(i,j,k)=1.0;
        }*/
        else
        {
            if(V_w_celltotal+V_a_celltotal > p->DXN[IP]*p->DYN[JP]*p->DZN[KP])
                cout<<"overflow"<<endl;
                
            a->vof(i,j,k)=V_w_celltotal/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            
            if(a->vof(i,j,k)>1.0)
            {
                cout<<"water overflow"<<endl;
                a->vof(i,j,k)=1.0;
            }
        }
    }
}