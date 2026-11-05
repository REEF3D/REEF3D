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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"interface_width.h"
#include"fluid_update_rheology.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"rheology_f.h"
#include"heaviside.h"

fluid_update_rheology::fluid_update_rheology(lexer *p) : ro1(p->W1), ro2(p->W3), visc2(p->W4)
{
    iter=0;
    iocheck = true;
    
    prheo = new rheology_f(p);

}

fluid_update_rheology::~fluid_update_rheology()
{
    delete prheo;
}

void fluid_update_rheology::start(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    const int gcval_ro = 1;
    const int gcval_visc = 1;

    double H_phi=0.0;
    p->volume1=0.0;
    p->volume2=0.0;
    
    if(p->count>iter)
    iocheck = 0;
    
    iter=p->count;

    // density, viscosity & volumes
    LOOP
    {  
        H_phi = heaviside(a->phi(i,j,k),interface_width(p,a->phi,i,j,k));
        

        a->ro(i,j,k) = ro1*H_phi + ro2*(1.0-H_phi);

        visc1 = prheo->viscosity(p,a,pgc,u,v,w);
        a->visc(i,j,k) = visc1*H_phi + visc2*(1.0-H_phi);

        if(p->flagsf4[IJK]>0)
        {
            p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H_phi*a->porosity(i,j,k);
            p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H_phi)*a->porosity(i,j,k);
        }
    }

    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    p->volume1 = pgc->globalsum(p->volume1);
    p->volume2 = pgc->globalsum(p->volume2);

    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
        cout<<"Volume 1: "<<p->volume1<<endl;
        cout<<"Volume 2: "<<p->volume2<<endl;
    }
    
    ++iocheck;
}
