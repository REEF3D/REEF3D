/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"net_void.h"


void net_void::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc
)
{
	
}

void net_void::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
}

void net_void::netForces
(
	double& Xme, double& Yme, double& Zme
)
{
	Xme = 0.0;
	Yme = 0.0;
	Zme = 0.0;
}