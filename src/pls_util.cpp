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
Author: Hans Bihs
--------------------------------------------------------------------*/
#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double particle_pls::hside(fdm* a)
{
    phival=fabs(a->phi(i,j,k));

        if(phival>epsi)
		Hval=0.0;

		if(phival<=epsi)
		Hval=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));

		return Hval;
}


void particle_pls::dgc_update(lexer* p,fdm* a,ghostcell* pgc)
{

    pgc->start1(p,a->u,14);
	pgc->start2(p,a->v,15);
	pgc->start3(p,a->w,16);
}

void particle_pls::vel_setback(lexer* p,fdm* a,ghostcell* pgc)
{/*
    pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);*/
}


