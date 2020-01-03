/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_boussinesq_void.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"slice1.h"
#include"slice2.h"
 
sflow_boussinesq_void::sflow_boussinesq_void(lexer* p, fdm2D *b) 
{
}

sflow_boussinesq_void::~sflow_boussinesq_void()
{
}

void sflow_boussinesq_void::ini(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    
}

void sflow_boussinesq_void::psi1(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{


}

void sflow_boussinesq_void::psi2(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    
}

