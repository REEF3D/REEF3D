/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"momentum_IM1.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"discrete.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"solver.h"
#include"ioflow.h"

void momentum_IM1::usource(lexer* p, fdm* a)
{
    count=0;
    ULOOP
    {
        a->M.p[count]+= 1.0/(UDT*CPOR1);

        a->F(i,j,k)+= a->M.p[count]*un(i,j,k)*(1.0/p->N54-1.0);

        a->M.p[count]/= p->N54;


    a->rhsvec.V[count] += (a->F(i,j,k) + a->gi)*PORVAL1 + un(i,j,k)/(UDT*CPOR1);

	++count;
    }
}

void momentum_IM1::vsource(lexer* p, fdm* a)
{
    count=0;
    VLOOP
    {

        a->M.p[count]+= 1.0/(VDT*CPOR2);

        a->G(i,j,k)+= a->M.p[count]*vn(i,j,k)*(1.0/p->N54-1.0);

        a->M.p[count]/= p->N54;


    a->rhsvec.V[count] += (a->G(i,j,k) + a->gj)*PORVAL2 + vn(i,j,k)/(VDT*CPOR2);

	++count;
    }
}

void momentum_IM1::wsource(lexer* p, fdm* a)
{
    count=0;
    WLOOP
    {

        a->M.p[count]+= 1.0/(WDT*CPOR3);

        a->H(i,j,k)+= a->M.p[count]*wn(i,j,k)*(1.0/p->N54-1.0);

        a->M.p[count]/= p->N54;

    a->rhsvec.V[count] +=(a->H(i,j,k) + a->gk)*PORVAL3 + wn(i,j,k)/(WDT*CPOR3);

	++count;
    }
}

void momentum_IM1::utimesave(lexer *p, fdm* a, ghostcell *pgc)
{
    ULOOP
    un(i,j,k)=a->u(i,j,k);
}

void momentum_IM1::vtimesave(lexer *p, fdm* a, ghostcell *pgc)
{
    VLOOP
    vn(i,j,k)=a->v(i,j,k);
}

void momentum_IM1::wtimesave(lexer *p, fdm* a, ghostcell *pgc)
{
    WLOOP
    wn(i,j,k)=a->w(i,j,k);
}

void momentum_IM1::clearrhs(lexer* p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	++count;
    }
}

void momentum_IM1::fillaij1(lexer *p, fdm *a, ghostcell *pgc, solver *psolv)
{
	pdisc->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
    usource(p,a);
}

void momentum_IM1::fillaij2(lexer *p, fdm *a, ghostcell *pgc, solver *psolv)
{
	pdisc->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
    vsource(p,a);
}

void momentum_IM1::fillaij3(lexer *p, fdm *a, ghostcell *pgc, solver *psolv)
{
	pdisc->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
    wsource(p,a);
}
