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

#include"rheology_v.h"

rheology_v::rheology_v() 
{

}

rheology_v::~rheology_v()
{
}

double rheology_v::viscosity(lexer*, fdm*, ghostcell*, field &u, field &v, field &w)
{
    return 0.0;
}

void rheology_v::u_source(lexer*, fdm*)
{
    
}

void rheology_v::v_source(lexer*, fdm*)
{
    
}

void rheology_v::w_source(lexer*, fdm*)
{
    
}

void rheology_v::filltau(lexer*, fdm*, ghostcell*)
{
}
