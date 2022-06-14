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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"LES_WALE.h"
#include"LES_filter_v.h"
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

LES_WALE::LES_WALE(lexer* p, fdm* a) : LES(p,a)
{
	gcval_sgs=24;
	c_wale=0.2;
    
    if(p->T21==0)
    pfilter = new LES_filter_v(p,a);
    
    if(p->T21==1)
    pfilter = new LES_filter_f1(p,a);
    
    if(p->T21==2)
    pfilter = new LES_filter_f2(p,a);
}

LES_WALE::~LES_WALE()
{
}

void LES_WALE::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow, vrans* pvrans)
{
    pfilter->start(p,a,pgc);
    
    LOOP
    a->eddyv(i,j,k) = pow(p->DXM*c_wale,2.0) *  (pow(magSqrSd(p,a), 3.0/2.0) / (pow(strainterm(p,a), 5.0) + pow(magSqrSd(p,a), 5.0/4.0)));

    pgc->start4(p,a->eddyv,gcval_sgs);
}

void LES_WALE::ktimesave(lexer* p, fdm* a, ghostcell *pgc)
{
}

void LES_WALE::etimesave(lexer* p, fdm* a, ghostcell *pgc)
{
}


