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

#include"wind_v.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"

wind_v::wind_v(lexer *p) 
{

}

wind_v::~wind_v()
{
}

void wind_v::wind_forcing_ini(lexer *p, ghostcell *pgc)
{
    
}

void wind_v::wind_forcing_nhf_x(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *F, slice &WL, slice &eta)
{
    
}

void wind_v::wind_forcing_nhf_y(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *G, slice &WL, slice &eta)
{
    
}