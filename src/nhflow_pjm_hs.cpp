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

#include"nhflow_pjm_hs.h"
#include"lexer.h"
#include"fdm_nhf.h" 
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"patchBC_interface.h"

#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)
#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0e20)
 
nhflow_pjm_hs::nhflow_pjm_hs(lexer* p, fdm_nhf *d, patchBC_interface *ppBC) : nhflow_gradient(p)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    gcval_press=540;  
}

nhflow_pjm_hs::~nhflow_pjm_hs()
{
}

void nhflow_pjm_hs::start(lexer*p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL, double *U, double *V, double *W, double alpha)
{
}

void nhflow_pjm_hs::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *U, double *P, double alpha)
{	
}

void nhflow_pjm_hs::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *V, double *P, double alpha)
{	 
}

void nhflow_pjm_hs::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *W, double *P, double alpha)
{
}
 
void nhflow_pjm_hs::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
}

void nhflow_pjm_hs::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
}

void nhflow_pjm_hs::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRY
    d->F[IJK] += PORVALNH*0.5*(d->ETAs(i,j)+d->ETAn(i-1,j))*fabs(p->W22)*
                (d->dfx(i,j) - d->dfx(i-1,j))/(p->DXN[IP]);
}

void nhflow_pjm_hs::vpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRY
	d->G[IJK] += PORVALNH*0.5*(d->ETAe(i,j)+d->ETAw(i,j-1))*fabs(p->W22)*
                 (d->dfy(i,j) - d->dfy(i,j-1))/(p->DYN[JP]);
}

void nhflow_pjm_hs::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
}
