/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"pjm_sig_hs.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
 
pjm_sig_hs::pjm_sig_hs(lexer* p, fdm *a, heat *&pheat, concentration *&pconc)
{
    if((p->F80==0||p->A10==5) && p->H10==0 && p->W30==0)
	pd = new density_f(p);
	
	if(p->F80==0 && p->H10==0 && p->W30==1)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0)
	pd = new density_conc(p,pconc);

    gcval_press=540;  

	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
}

pjm_sig_hs::~pjm_sig_hs()
{
}

void pjm_sig_hs::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
}

void pjm_sig_hs::ucorr(lexer* p, fdm* a, field& uvel,double alpha)
{	
}

void pjm_sig_hs::vcorr(lexer* p, fdm* a, field& vvel,double alpha)
{	 
}

void pjm_sig_hs::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{
}
 
void pjm_sig_hs::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
}
 
void pjm_sig_hs::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
}

void pjm_sig_hs::upgrad(lexer*p,fdm* a)
{
    ULOOP
	a->F(i,j,k) -= PORVAL1*fabs(p->W22)*(a->eta(i+1,j)-a->eta(i,j))/p->DXP[IP];
}

void pjm_sig_hs::vpgrad(lexer*p,fdm* a)
{
    VLOOP
	a->G(i,j,k) -= PORVAL2*fabs(p->W22)*(a->eta(i,j+1)-a->eta(i,j))/p->DYP[JP];
}

void pjm_sig_hs::wpgrad(lexer*p,fdm* a)
{
    WLOOP
	a->H(i,j,k) -= a->gk*PORVAL3;
}





