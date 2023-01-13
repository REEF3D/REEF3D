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

#include"vrans_veg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_veg::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
	int qn;
    double zmin,zmax,slope;
    double xs,xe;
	
	
	ALOOP
	{
	N(i,j,k)=0.0;
	D(i,j,k)=0.0;
	Cd(i,j,k)=0.0;
    a->porosity(i,j,k)=1.0;
	}
	
	pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,N,1);
	pgc->start4a(p,D,1);
	pgc->start4a(p,Cd,1);
	
	
	
	// Box
    for(qn=0;qn<p->B310;++qn)
    ALOOP
	if(p->XN[IP]>=p->B310_xs[qn] && p->XN[IP]<p->B310_xe[qn]
	&& p->YN[JP]>=p->B310_ys[qn] && p->YN[JP]<p->B310_ye[qn]
	&& p->ZN[KP]>=p->B310_zs[qn] && p->ZN[KP]<p->B310_ze[qn])
	{
	N(i,j,k) = p->B310_N[qn];
	D(i,j,k) = p->B310_D[qn];
	Cd(i,j,k) = p->B310_Cd[qn];
    
    if(p->B308==1)
    a->porosity(i,j,k) =  1.0 - (p->B310_N[qn]*PI*pow(p->B310_D[qn],2.0)*0.25); 
	}
    
    if(p->mpirank==0)
    for(qn=0;qn<p->B310;++qn)
    cout<<"VEG :  "<<(p->B310_N[qn]*PI*pow(p->B310_D[qn],2.0)*0.25)<<endl;

    
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,N,1);
	pgc->start4a(p,D,1);
	pgc->start4a(p,Cd,1);
}

