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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"

void  ioflow_f::isource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
    double porousterm;

	// Darcy Porosity
	count=0;
    if(p->B240>0 && p->B241==1)
    ULOOP
	{
		
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->u(i,j,k) + 0.5*p->B240_C[n]*a->u(i,j,k)*fabs(a->u(i,j,k));
		}
	
    a->rhsvec.V[count] -= porousterm;
	++count;
	}
	
	//VRANS
    pvrans->u_source(p,a);
}

void  ioflow_f::jsource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
    double porousterm;

	count=0;
    if(p->B240>0 && p->B242==1)
    VLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->v(i,j,k) + 0.5*p->B240_C[n]*a->v(i,j,k)*fabs(a->v(i,j,k));
		}
	
    a->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
    pvrans->v_source(p,a);
}

void  ioflow_f::ksource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
    double porousterm;
	
	count=0;
    if(p->B240>0 && p->B243==1)
    WLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->w(i,j,k) + 0.5*p->B240_C[n]*a->w(i,j,k)*fabs(a->w(i,j,k));
		}

    a->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
    pvrans->w_source(p,a);
}


