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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sflow_turb_parabolic.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

#define HPIJ (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

sflow_turb_parabolic::sflow_turb_parabolic(lexer* p) : sflow_turb_io_void(p)
{
}

sflow_turb_parabolic::~sflow_turb_parabolic()
{
}

void sflow_turb_parabolic::start(lexer *p, fdm2D *b, ghostcell *pgc, sflow_convection *pconvec, sflow_diffusion *pdiff, solver2D *psolv, ioflow *pflow)
{
    double dudx,dvdy,dudy,dvdx;
    double alpha_t,Ustar;
    double manning,cf;
    
    alpha_t = p->A262;
    
	SLICELOOP4
    {
    manning = pow(b->ks(i,j),1.0/6.0)/26.0;
    
    cf = pow(manning,2.0)*9.81/pow(HPIJ,1.0/3.0);
    
    Ustar = sqrt(cf*(b->P(i,j)*b->P(i,j) + b->Q(i,j)*b->Q(i,j)));
    
    b->eddyv(i,j) = alpha_t*Ustar*b->hp(i,j);
    }
    
    pgc->gcsl_start4(p,b->eddyv,24);
}

void sflow_turb_parabolic::ktimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{

}

void sflow_turb_parabolic::etimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{

}
