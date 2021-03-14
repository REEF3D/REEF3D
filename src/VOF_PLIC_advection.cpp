/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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


void VOF_PLIC::advectPlane
(
	fdm* a, 
	lexer* p, 
	const double Q1, 
	const double Q2, 
	int sweep
)
{
	//- According to Wang (24,25)
	
	if (sweep == 0)
	{
		nx(i, j, k) /= (1.0 + (Q2 - Q1));
		alpha(i, j, k) += nx(i, j, k)*Q1;
	}
	else if (sweep == 1)
	{
		ny(i, j, k) /= (1.0 + (Q2 - Q1));
		alpha(i, j, k) += ny(i, j, k)*Q1;
	}
	else
	{
		nz(i, j, k) /= (1.0 + (Q2 - Q1));
		alpha(i, j, k) += nz(i, j, k)*Q1;
	}
}


void VOF_PLIC::calcFlux(fdm* a, lexer* p, double& Q1, double& Q2, int sweep)
{
	double vell, velr;
	
	if (sweep == 0)
	{
		vell = a->u(i-1, j, k);
		velr = a->u(i, j, k);

		Q1 = ((velr - vell)*vell*p->dt*p->dt/(2.0*p->DXN[IP]) + vell*p->dt)/p->DXN[IP];
		Q2 = ((velr - vell)*velr*p->dt*p->dt/(2.0*p->DXN[IP]) + velr*p->dt)/p->DXN[IP];			
	}
	else if (sweep == 1)
	{
		vell = a->v(i, j-1, k);
		velr = a->v(i, j, k);
		
		Q1 = ((velr - vell)*vell*p->dt*p->dt/(2.0*p->DYN[JP]) + vell*p->dt)/p->DYN[JP];
		Q2 = ((velr - vell)*velr*p->dt*p->dt/(2.0*p->DYN[JP]) + velr*p->dt)/p->DYN[JP];	
	}
	else
	{
		vell = a->w(i, j, k-1);
		velr = a->w(i, j, k);

		Q1 = ((velr - vell)*vell*p->dt*p->dt/(2.0*p->DZN[KP]) + vell*p->dt)/p->DZN[KP];
		Q2 = ((velr - vell)*velr*p->dt*p->dt/(2.0*p->DZN[KP]) + velr*p->dt)/p->DZN[KP];
	}
}
	
