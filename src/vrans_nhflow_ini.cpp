/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"vrans_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void vrans_nhflow::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
    // ************************
    // start raycast
    geometry_ini(p, d, pgc);
    objects_create_forcing(p, pgc);
    ray_cast(p, d, pgc, d->PORSTRUC);
    reini_RK2(p, d, pgc, d->PORSTRUC);
    // ************************
    
    /*
	int qn;
    double zmin,zmax,slope;
    double xs,xe,ys,ye;
	
	LOOP
	{
	a->porosity(i,j,k)=1.0;
	a->porpart(i,j,k)=0.01;
	alpha(i,j,k)=0.0;
	beta(i,j,k)=0.0;
	}
	
	pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,a->porpart,1);
	pgc->start4a(p,alpha,1);
	pgc->start4a(p,beta,1);
	
	
	// Box
    for(qn=0;qn<p->B270;++qn)
    ALOOP
	if(p->XN[IP]>=p->B270_xs[qn] && p->XN[IP]<p->B270_xe[qn]
	&& p->YN[JP]>=p->B270_ys[qn] && p->YN[JP]<p->B270_ye[qn]
	&& p->ZN[KP]>=p->B270_zs[qn] && p->ZN[KP]<p->B270_ze[qn])
	{
	a->porosity(i,j,k)= p->B270_n[qn];
	a->porpart(i,j,k) = p->B270_d50[qn];
	alpha(i,j,k) = p->B270_alpha[qn];
	beta(i,j,k) = p->B270_beta[qn];
	}
    
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,a->porpart,1);
	pgc->start4a(p,alpha,1);
	pgc->start4a(p,beta,1);
    

    */
}

