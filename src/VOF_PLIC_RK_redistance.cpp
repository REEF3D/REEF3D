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
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"lexer.h"
#include"gradient.h"

void VOF_PLIC::RK_redistance
(
    fdm* a,
    lexer* p,
    ghostcell* pgc
)
{
    LOOP
    {
        if(a->vof(i,j,k)>0.001 && a->vof(i,j,k)<0.999)
            reconstructPlane_alt(a,p,a->vof);
    }
    else
    {
        nx(i,j,k)=2.0;
        ny(i,j,k)=0.0;
        nz(i,j,k)=2.0;
        alpha(i,j,k)=1E20;
    }
    phiaux(i,j,k)=1E05;
    
    pgc->start4(p,nx,1);
    pgc->start4(p,ny,1);
    pgc->start4(p,nz,1);
    pgc->start4(p,alpha,1);
    pgc->start4(p,phiaux,1);

    LOOP
    {
        if(a->vof(i,j,k)>0.001 && a->vof(i,j,k)<0.999)
            redistancePhiByPlane_Bonn(a,p);
    }
    pgc->start4(p,phiaux,1);
    LOOP
    {
        if(fabs(phiaux(i,j,k))<1E04)
        {
            a->phi(i,j,k)=phiaux(i,j,k);
        }
        else
        {
            if(a->vof(i,j,k)>=0.5)
                a->phi(i,j,k)=0.5;
            else
                a->phi(i,j,k)=-0.5;
            
               
        }
    }
    pgc->start4(p,a->phi,1);

}