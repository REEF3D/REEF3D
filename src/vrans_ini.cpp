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

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_f::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
	int qn;
    double zmin,zmax,slope;
    double xs,xe,ys,ye;
	
	ALOOP
	{
        a->porosity(i,j,k)=1.0;
        porpart(i,j,k)=0.01;
        alpha(i,j,k)=0.0;
        beta(i,j,k)=0.0;
	}
	
	pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,porpart,1);
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
                porpart(i,j,k) = p->B270_d50[qn];
                alpha(i,j,k) = p->B270_alpha[qn];
                beta(i,j,k) = p->B270_beta[qn];
            }
    
    // Vertical Cylinder
    for(qn=0;qn<p->B274;++qn)
        ALOOP
        {
            double  r = sqrt( pow(p->XP[IP]-p->B274_xc[qn],2.0)+pow(p->YP[JP]-p->B274_yc[qn],2.0));
            
            if(r<=p->B274_r[qn] && p->pos_z()>p->B274_zs[qn] && p->pos_z()<=p->B274_ze[qn])
            {
                a->porosity(i,j,k)= p->B274_n[qn];
                porpart(i,j,k) = p->B274_d50[qn];
                alpha(i,j,k) = p->B274_alpha[qn];
                beta(i,j,k) = p->B274_beta[qn];
            }
        }

	
	// Wedge x-dir
    for(qn=0;qn<p->B281;++qn)
    {
		zmin=MIN(p->B281_zs[qn],p->B281_ze[qn]);
        
        if(p->B281_xs[qn]<=p->B281_xe[qn])
        {
            xs = p->B281_xs[qn];
            xe = p->B281_xe[qn];
        }
            
        if(p->B281_xs[qn]>p->B281_xe[qn])
        {
            xs = p->B281_xe[qn];
            xe = p->B281_xs[qn];
        }

		slope=(p->B281_ze[qn]-p->B281_zs[qn])/(p->B281_xe[qn]-p->B281_xs[qn]);

		ALOOP
            if(p->pos_x()>=xs && p->pos_x()<xe
            && p->pos_y()>=p->B281_ys[qn] && p->pos_y()<p->B281_ye[qn]
            && p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_x()-p->B281_xs[qn])+p->B281_zs[qn] )
            {
                a->porosity(i,j,k)=p->B281_n[qn];
                porpart(i,j,k) =p->B281_d50[qn];
                alpha(i,j,k) = p->B281_alpha[qn];
                beta(i,j,k) = p->B281_beta[qn];
            }
    }
    
    // Wedge y-dir
    for(qn=0;qn<p->B282;++qn)
    {
		zmin=MIN(p->B282_zs[qn],p->B282_ze[qn]);
        
        if(p->B282_xs[qn]<=p->B282_xe[qn])
        {
            ys = p->B282_ys[qn];
            ye = p->B282_ye[qn];
        }
            
        if(p->B282_xs[qn]>p->B282_xe[qn])
        {
            ys = p->B282_ye[qn];
            ye = p->B282_ys[qn];
        }

		slope=(p->B282_ze[qn]-p->B282_zs[qn])/(p->B282_ye[qn]-p->B282_ys[qn]);

		ALOOP
            if(p->pos_x()>=p->B282_xs[qn] && p->pos_x()<p->B282_xe[qn]
            && p->pos_y()>=ys && p->pos_y()<ye
            && p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_y()-p->B282_ys[qn])+p->B282_zs[qn] )
            {
                a->porosity(i,j,k)=p->B282_n[qn];
                porpart(i,j,k) =p->B282_d50[qn];
                alpha(i,j,k) = p->B282_alpha[qn];
                beta(i,j,k) = p->B282_beta[qn];
            }
    }
    
    // Plate x-dir
    for(qn=0;qn<p->B291;++qn)
    {
		zmin=MIN(p->B291_zs[qn],p->B291_ze[qn]);
        zmin=MAX(p->B291_zs[qn],p->B291_ze[qn]);
        
        if(p->B291_xs[qn]<=p->B291_xe[qn])
        {
            xs = p->B291_xs[qn];
            xe = p->B291_xe[qn];
        }
            
        if(p->B291_xs[qn]>p->B291_xe[qn])
        {
            xs = p->B291_xe[qn];
            xe = p->B291_xs[qn];
        }

		slope=(p->B291_ze[qn]-p->B291_zs[qn])/(p->B291_xe[qn]-p->B291_xs[qn]);

		ALOOP
            if(p->pos_x()>=xs && p->pos_x()<xe
            && p->pos_y()>=p->B291_ys[qn] && p->pos_y()<p->B291_ye[qn]
            
            && p->pos_z()>=zmin 
            && p->pos_z()<=zmax 
            
            && p->pos_z()<slope*(p->pos_x()-p->B291_xs[qn])+p->B291_zs[qn]+p->B291_d[qn] //lower
            && p->pos_z()>slope*(p->pos_x()-p->B291_xs[qn])+p->B291_zs[qn]) // upper
            {
                a->porosity(i,j,k)=p->B291_n[qn];
                porpart(i,j,k) =p->B291_d50[qn];
                alpha(i,j,k) = p->B291_alpha[qn];
                beta(i,j,k) = p->B291_beta[qn];
            }
    }
    
    
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,porpart,1);
	pgc->start4a(p,alpha,1);
	pgc->start4a(p,beta,1);
    
    
    // Sediment
    if(p->S10==2)
        sed_update(p,a,pgc);
}

