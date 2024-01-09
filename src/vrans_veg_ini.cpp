/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
    double xs,xe,ys,ye;
	
	
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
    
    // Wedge x-dir
    for(qn=0;qn<p->B321;++qn)
    {
		zmin=MIN(p->B321_zs[qn],p->B321_ze[qn]);
        
            if(p->B321_xs[qn]<=p->B321_xe[qn])
            {
            xs = p->B321_xs[qn];
            xe = p->B321_xe[qn];
            }
            
            if(p->B321_xs[qn]>p->B321_xe[qn])
            {
            xs = p->B321_xe[qn];
            xe = p->B321_xs[qn];
            }

		slope=(p->B321_ze[qn]-p->B321_zs[qn])/(p->B321_xe[qn]-p->B321_xs[qn]);

		ALOOP
		if(p->pos_x()>=xs && p->pos_x()<xe
		&& p->pos_y()>=p->B321_ys[qn] && p->pos_y()<p->B321_ye[qn]
		&& p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_x()-p->B321_xs[qn])+p->B321_zs[qn] )
		{
		N(i,j,k) = p->B321_N[qn];
         D(i,j,k) = p->B321_D[qn];
        Cd(i,j,k) = p->B321_Cd[qn];
		}
    }
    
    // Wedge y-dir
    for(qn=0;qn<p->B322;++qn)
    {
		zmin=MIN(p->B322_zs[qn],p->B322_ze[qn]);
        
            if(p->B322_xs[qn]<=p->B322_xe[qn])
            {
            ys = p->B322_ys[qn];
            ye = p->B322_ye[qn];
            }
            
            if(p->B322_xs[qn]>p->B322_xe[qn])
            {
            ys = p->B322_ye[qn];
            ye = p->B322_ys[qn];
            }

		slope=(p->B322_ze[qn]-p->B322_zs[qn])/(p->B322_ye[qn]-p->B322_ys[qn]);

		ALOOP
		if(p->pos_x()>=p->B322_xs[qn] && p->pos_x()<p->B322_xe[qn]
		&& p->pos_y()>=ys && p->pos_y()<ye
		&& p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_y()-p->B322_ys[qn])+p->B322_zs[qn] )
		{
		N(i,j,k) = p->B322_N[qn];
         D(i,j,k) = p->B322_D[qn];
        Cd(i,j,k) = p->B322_Cd[qn];
		}
    }
    

    
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,N,1);
	pgc->start4a(p,D,1);
	pgc->start4a(p,Cd,1);
}

