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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::interface(lexer *p, bool final)
{
	p->ufb=Us;
	p->vfb=Vs;
	p->wfb=Ws;
	
	p->pfb=Ps;
	p->qfb=Qs;
	p->rfb=Rs;
	
	if (final == true)
	{
		p->ufbn=p->ufbi;
		p->vfbn=p->vfbi;
		p->wfbn=p->wfbi;
		
		p->pfbn=p->pfbi;
		p->qfbn=p->qfbi;
		p->rfbn=p->rfbi;
	}
	
	p->ufbi=Ue+Uext;
	p->vfbi=Ve+Vext;
	p->wfbi=We+Wext;
	
	p->pfbi=Pe+Pext;
	p->qfbi=Qe+Qext;
	p->rfbi=Re+Rext;
	
	p->xg=xg;
	p->yg=yg;
	p->zg=zg;
	
	p->xgn=xgn;
	p->ygn=ygn;
	p->zgn=zgn;
	
	p->phi_fb=phi;
	p->theta_fb=theta;
	p->psi_fb=psi;
}

void sixdof_gc::maxvel(lexer *p, fdm *a, ghostcell *pgc)
{
	double uvel,vvel,wvel;
	
	p->ufbmax=p->vfbmax=p->wfbmax=0.0;
	uvel=vvel=wvel=0.0;
	
	ALOOP
	{
	uvel = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
	vvel = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
	wvel = p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi;
	
	p->ufbmax = MAX(p->ufbmax,fabs(uvel));
	p->vfbmax = MAX(p->vfbmax,fabs(vvel));
	p->wfbmax = MAX(p->wfbmax,fabs(wvel));
	}
	
	
}

