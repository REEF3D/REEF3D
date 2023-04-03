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
Authors: Tobias Martin, Ahmet Soydan, Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::solid_forcing_ini(lexer *p, fdm *a)
{
    // Initialise floating fields
     ULOOP
     a->fbh1(i,j,k) = Hsolidface(p,a,1,0,0);

     VLOOP
     a->fbh2(i,j,k) = Hsolidface(p,a,0,1,0);

     WLOOP
     a->fbh3(i,j,k) = Hsolidface(p,a,0,0,1);

     LOOP
     a->fbh4(i,j,k) = Hsolidface(p,a,0,0,0);

     start1(p,a->fbh1,10);
     start2(p,a->fbh2,11);
     start3(p,a->fbh3,12);
     start4(p,a->fbh4,40);
     
     
     // ghostcell update
    gcdf_update(p,a);
    
}