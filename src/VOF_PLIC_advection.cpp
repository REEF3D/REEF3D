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
Author: Tobias Martin
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
	int sweep
)
{
	//- According to Wang (24,25)
	double N_sum;
	if (sweep == 0)
	{
		nx(i, j, k) /= (1.0 + (Q2(i,j,k) - Q1(i,j,k))/p->DXN[IP]);
        N_sum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
        nx(i,j,k)/=N_sum;
		alpha(i, j, k) += nx(i, j, k)*Q1(i,j,k);
	}
	else if (sweep == 1)
	{
		ny(i, j, k) /= (1.0 + (Q2(i,j,k) - Q1(i,j,k))/p->DYN[JP]);
        N_sum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
        ny(i,j,k)/=N_sum;
		alpha(i, j, k) += ny(i, j, k)*Q1(i,j,k);
	}
	else
	{
		nz(i, j, k) /= (1.0 + (Q2(i,j,k) - Q1(i,j,k))/p->DZN[KP]);
        N_sum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
        nz(i,j,k)/=N_sum;
		alpha(i, j, k) += nz(i, j, k)*Q1(i,j,k);
	}
}


void VOF_PLIC::calcFlux(fdm* a, lexer* p, int sweep)
{
	double vell, velr;
	
	if (sweep == 0)
	{
		vell = a->u(i-1, j, k);
		velr = a->u(i, j, k);

		Q1(i,j,k) = (velr - vell)*vell*p->dt*p->dt/(2.0*p->DXN[IP]) + vell*p->dt;
		Q2(i,j,k) = (velr - vell)*velr*p->dt*p->dt/(2.0*p->DXN[IP]) + velr*p->dt;
	}
	else if (sweep == 1)
	{
		vell = a->v(i, j-1, k);
		velr = a->v(i, j, k);
		
		Q1(i,j,k) = (velr - vell)*vell*p->dt*p->dt/(2.0*p->DYN[JP]) + vell*p->dt;
		Q2(i,j,k) = (velr - vell)*velr*p->dt*p->dt/(2.0*p->DYN[JP]) + velr*p->dt;
	}
	else
	{
		vell = a->w(i, j, k-1);
		velr = a->w(i, j, k);

		Q1(i,j,k) = (velr - vell)*vell*p->dt*p->dt/(2.0*p->DZN[KP]) + vell*p->dt;
		Q2(i,j,k) = (velr - vell)*velr*p->dt*p->dt/(2.0*p->DZN[KP]) + velr*p->dt;
	}
}
	
