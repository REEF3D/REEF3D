/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"

void iowave::isource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;

    double porousterm;
    double ep=1.0e-10*p->DXM;
	count=0;
	if(p->B240>0 && p->B241==1)
    ULOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos1_x() > p->B240_xs[n]+ep && p->pos1_x() <= p->B240_xe[n]+ep)
			if(p->pos1_y() > p->B240_ys[n]+ep && p->pos1_y() <= p->B240_ye[n]+ep)
			if(p->pos1_z() > p->B240_zs[n]+ep && p->pos1_z() <= p->B240_ze[n]+ep)
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->u(i,j,k) + 0.5*p->B240_C[n]*a->u(i,j,k)*fabs(a->u(i,j,k));
		}

    a->rhsvec.V[count] -= porousterm;
	++count;
	}

	//VRANS
    pvrans->u_source(p,a);
}

void iowave::jsource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;

    double porousterm;
	double ep=1.0e-10*p->DXM;
	count=0;
	if(p->B240>0 && p->B242==1)
    VLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos2_x() > p->B240_xs[n]+ep && p->pos2_x() <= p->B240_xe[n]+ep)
			if(p->pos2_y() > p->B240_ys[n]+ep && p->pos2_y() <= p->B240_ye[n]+ep)
			if(p->pos2_z() > p->B240_zs[n]+ep && p->pos2_z() <= p->B240_ze[n]+ep)
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->v(i,j,k) + 0.5*p->B240_C[n]*a->v(i,j,k)*fabs(a->v(i,j,k));
		}

    a->rhsvec.V[count] -= porousterm;
	++count;
	}

    //VRANS
    pvrans->v_source(p,a);
}

void iowave::ksource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;

    double porousterm;


    double ep=1.0e-10*p->DXM;
	count=0;
	if(p->B240>0 && p->B243==1)
    WLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos3_x() > p->B240_xs[n]+ep && p->pos3_x() <= p->B240_xe[n]+ep)
			if(p->pos3_y() > p->B240_ys[n]+ep && p->pos3_y() <= p->B240_ye[n]+ep)
			if(p->pos3_z() > p->B240_zs[n]+ep && p->pos3_z() <= p->B240_ze[n]+ep)
            {
            //cout<<"k: "<<k<<" pos_z: "<<p->pos3_z()<<endl;
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->w(i,j,k) + 0.5*p->B240_C[n]*a->w(i,j,k)*fabs(a->w(i,j,k));
            }
		}

    a->rhsvec.V[count] -= porousterm;
	++count;
	}

    //VRANS
    pvrans->w_source(p,a);
}

void iowave::isource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
	SLICELOOP1
	b->F(i,j)=0.0;
}

void iowave::jsource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
	SLICELOOP2
	b->G(i,j)=0.0;
}
