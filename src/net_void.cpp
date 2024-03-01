/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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
	ghostcell *pgc,
    double alpha,
    Eigen::Matrix3d quatRotMat
)
{
	
}

void net_void::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
}

void net_void::netForces
(
    lexer *p,
	double& Xne, double& Yne, double& Zne,
	double& Kne, double& Mne, double& Nne
)
{
	Xne = 0.0;
	Yne = 0.0;
	Zne = 0.0;
	Kne = 0.0;
	Mne = 0.0;
	Nne = 0.0;
}
