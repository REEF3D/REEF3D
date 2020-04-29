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

#include"LES_smagorinsky.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"solver.h"
#include"diffusion.h"
#include"ioflow.h"
#include"convection.h"

LES_smagorinsky::LES_smagorinsky(lexer* p, fdm* a) : LES(p,a)
{
	gcval_sgs=24;
	c_sgs=0.2;
}

LES_smagorinsky::~LES_smagorinsky()
{
}

void LES_smagorinsky::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    LOOP
    a->eddyv(i,j,k) = pow(p->dx*c_sgs,2.0) * strainterm(p,a);
    
    LOOP
    a->visctot(i,j,k) = a->visc(i,j,k) + a->eddyv(i,j,k);

    pgc->start4(p,a->eddyv,gcval_sgs);
    pgc->start4(p,a->visctot,1);
}

void LES_smagorinsky::ktimesave(lexer* p, fdm* a, ghostcell *pgc)
{
}

void LES_smagorinsky::etimesave(lexer* p, fdm* a, ghostcell *pgc)
{
}


