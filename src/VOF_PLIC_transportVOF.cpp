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

void VOF_PLIC::vof_transport_COSMIC2D
(
    fdm* a,
    lexer* p,
    int nSweep,
    int sweep
)
{
    if(nSweep<2 && sweep==0)
    {   
        LOOP
        {
            Flux_x(i,j,k)=(Vn_p(i,j,k)-Vn_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

            F_x(i,j,k)=F_n(i,j,k)
                        -Flux_x(i,j,k)
                        +p->dt*F_n(i,j,k)*(a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP];
                    
            if(nSweep==1)
            {
                Crossflux_zx(i,j,k)=(Vz_p(i,j,k)-Vz_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
        }
    }
    
    if(nSweep<2 && sweep==2)
    {
        LOOP
        {
            Flux_z(i,j,k)=(Vn_p(i,j,k)-Vn_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
        
            F_z(i,j,k)=F_n(i,j,k)
                        -Flux_z(i,j,k)
                        +p->dt*F_n(i,j,k)*(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP];
                    
            if(nSweep==1)
            {
                Crossflux_xz(i,j,k)=(Vx_p(i,j,k)-Vx_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
        }
    }
    
    if(nSweep==2)
    {
        LOOP
        {
            if(sweep==0)
            {
                Crossflux_xz(i,j,k)=(Vx_p(i,j,k)-Vx_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
            else if(sweep==2)
            {
                Crossflux_zx(i,j,k)=(Vz_p(i,j,k)-Vz_m(i,j,k))/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            }
        
            F_new(i,j,k)=F_n(i,j,k)-0.5*(Flux_x(i,j,k)+Crossflux_zx(i,j,k)+Flux_z(i,j,k)+Crossflux_xz(i,j,k));
        }               
    }
    
}