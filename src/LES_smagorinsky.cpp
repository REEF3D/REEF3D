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

#include"LES_smagorinsky.h"
#include"LES_filter_box.h"
#include"LES_filter_f1.h"
#include"LES_filter_f2.h"
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
    gcval_u1=10;
	gcval_v1=11;
	gcval_w1=12;
    
	gcval_sgs=24;
	c_sgs=0.2;
    
    if(p->T21==0)
    pfilter = new LES_filter_box(p,a);
    
    if(p->T21==1)
    pfilter = new LES_filter_f1(p,a);
    
    if(p->T21==2)
    pfilter = new LES_filter_f2(p,a);
}

LES_smagorinsky::~LES_smagorinsky()
{
}

void LES_smagorinsky::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow, vrans* pvrans)
{
    pfilter->start(p,a,pgc,uprime,vprime,wprime,gcval_u1);
    pfilter->start(p,a,pgc,uprime,vprime,wprime,gcval_v1);
    pfilter->start(p,a,pgc,uprime,vprime,wprime,gcval_w1);
	
    
    LOOP
    a->eddyv(i,j,k) = pow(c_sgs,2.0) * pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP],2.0/3.0) * sqrt(2.0) * strainterm(p,uprime,vprime,wprime);

//		a->eddyv(i,j,k) = pow(p->DXM*c_sgs,2.0) * sqrt(2.0) * strainterm(p,uprime,vprime,wprime);
//    a->eddyv(i,j,k) = pow(p->DXM*c_sgs,2.0) * sqrt(2.0) * strainterm(p,a);

    pgc->start4(p,a->eddyv,gcval_sgs);
}

void LES_smagorinsky::ktimesave(lexer* p, fdm* a, ghostcell *pgc)
{
}

void LES_smagorinsky::etimesave(lexer* p, fdm* a, ghostcell *pgc)
{
}


