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
#include"sflow_turb_prandtl.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_turb_prandtl::sflow_turb_prandtl(lexer* p) : sflow_turb_io_void(p)
{
}

sflow_turb_prandtl::~sflow_turb_prandtl()
{
}

void sflow_turb_prandtl::start(lexer *p, fdm2D *b, ghostcell *pgc, sflow_convection *pconvec, sflow_diffusion *pdiff, solver2D *psolv, ioflow *pflow)
{
    double dudx,dvdy,dudy,dvdx;
    double refl;
    
	SLICELOOP4
    {
    dudx=dvdy=dudy=dvdx=0.0;
    
    if(p->flagslice1[IJ]>0 && p->flagslice1[Im1J]>0)
    dudx = (b->P(i,j) - b->P(i-1,j))/(p->DXM);
    
    if(p->flagslice2[IJ]>0 && p->flagslice2[IJm1]>0)
    dvdy = (b->Q(i,j) - b->Q(i,j-1))/(p->DXM);
    
    if(p->flagslice1[IJp1]>0 && p->flagslice1[Im1Jp1]>0 && p->flagslice1[IJm1]>0 && p->flagslice1[Im1Jm1]>0)
    dudy = (0.5*(b->P(i,j+1)+b->P(i-1,j+1)) - 0.5*(b->P(i,j-1)+b->P(i-1,j-1)))/(2.0*p->DXM);
    
    if(p->flagslice2[Ip1J]>0 && p->flagslice2[Ip1Jm1]>0 && p->flagslice2[Im1J]>0 && p->flagslice2[Im1Jm1]>0)
    dvdx = (0.5*(b->Q(i+1,j)+b->Q(i+1,j-1)) - 0.5*(b->Q(i-1,j)+b->Q(i-1,j-1)))/(2.0*p->DXM);
    
    refl = 0.4*p->A261*b->hp(i,j);
 
    b->eddyv(i,j) = pow(refl,2.0)*sqrt(2.0*pow(dudx,2.0) + 2.0*pow(dvdy,2.0) + pow(dudy+dvdx,2.0));
    }
    
    pgc->gcsl_start4(p,b->eddyv,24);
}

void sflow_turb_prandtl::ktimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{

}

void sflow_turb_prandtl::etimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{

}
