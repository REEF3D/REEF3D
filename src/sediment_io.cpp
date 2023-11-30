/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/
#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"bedshear.h"
#include"suspended.h"

double sediment_f::bedshear_point(lexer *p, fdm *a,ghostcell *pgc)
{
	double tau_eff = s->tau_eff(i,j);
    
	return tau_eff;
}

void sediment_f::fill_bedk(lexer *p, fdm *a,ghostcell *pgc)
{
    SLICELOOP4
    s->bedk(i,j)=0;
    
    SLICELOOP4
    KLOOP
    PBASECHECK
    if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
    s->bedk(i,j)=k+1;
}

void sediment_f::fill_PQ_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    double zval,xip,yip;
    
    SLICELOOP1
    {
    k=s->bedk(i,j);
    
    xip= p->XN[IP1];
	yip= p->YP[JP];
    zval = 0.5*(s->bedzh(i,j)+s->bedzh(i+1,j)) + 1.6*p->DZN[k];
    
    s->P(i,j) = a->P(i,j) = p->ccipol1_a(a->u,xip,yip,zval);
    }
    
    SLICELOOP2
    {
    k=s->bedk(i,j);
    
    xip= p->XP[IP];
	yip= p->YN[JP1];
    zval = 0.5*(s->bedzh(i,j)+s->bedzh(i,j+1)) + 1.6*p->DZN[k];
    
    s->Q(i,j) = a->Q(i,j)  = p->ccipol2_a(a->v,xip,yip,zval);
    }
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);
    
    pgc->gcsl_start1(p,a->P,10);
	pgc->gcsl_start2(p,a->Q,11);
}

void sediment_f::fill_PQ_sflow(lexer *p, fdm2D *b,ghostcell *pgc,slice &P, slice &Q)
{
    SLICELOOP1
    s->P(i,j) = P(i,j);
    
    SLICELOOP2
    s->Q(i,j) = Q(i,j);
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);
}

double sediment_f::qbeval(int ii, int jj)
{
    double val;

    val=s->qbe(ii,jj);

    return val;
}

void sediment_f::qbeget(int ii, int jj, double val)
{
    s->qbe(ii,jj)=val;
}

double sediment_f::bedzhval(int ii, int jj)
{
    double val;

    val=s->bedzh(ii,jj);

    return val;
}

void sediment_f::ctimesave(lexer *p, fdm* a)
{
    psusp->ctimesave(p,a);
}
