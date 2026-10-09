/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_f::initialize_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    structures(p,a);
    
    // continuum sediment bed (S 10 2 with the Exner sediment module, Q 10 0)
    if(p->S10==2 && p->Q10==0)
    sediment_bed(p,a);
    
    exchange(p,a,pgc);
}

// cell values of a porous zone; the Darcy-Forchheimer coefficients (van Gent 1995, Jensen et al. 2014)
// A = porA nu,  porA = alpha (1-n)^2/(n^3 d50^2)
// B = porB,     porB = beta (1 + 7.5/KC) (1-n)/(n^3 d50)
void vrans_f::set_cell(lexer *p, fdm *a, double n, double d50, double alpha, double beta)
{
    a->porosity(i,j,k) = n;
    a->porpart(i,j,k) = d50;
    a->porA(i,j,k) = 0.0;
    a->porB(i,j,k) = 0.0;
    
    if(n<1.0)
    {
    a->porA(i,j,k) = alpha*pow(1.0-n,2.0)/(pow(n,3.0)*d50*d50);
    a->porB(i,j,k) = beta*(1.0 + 7.5/Cval)*(1.0-n)/(pow(n,3.0)*d50);
    }
}

void vrans_f::exchange(lexer *p, fdm *a, ghostcell *pgc)
{
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,a->porpart,1);
	pgc->start4a(p,a->porA,1);
	pgc->start4a(p,a->porB,1);
}

// reset to open flow, then the porous structures B 270 - B 291
// (called before every sediment update, so the structures are kept)
void vrans_f::structures(lexer *p, fdm *a)
{
	int qn;
    double zmin,zmax,slope;
    double xs,xe,ys,ye;
	
	BASELOOP
	set_cell(p,a,1.0,0.01,0.0,0.0);
	
	// Box
    for(qn=0;qn<p->B270;++qn)
    LOOP
	if(p->XN[IP]>=p->B270_xs[qn] && p->XN[IP]<p->B270_xe[qn]
	&& p->YN[JP]>=p->B270_ys[qn] && p->YN[JP]<p->B270_ye[qn]
	&& p->ZN[KP]>=p->B270_zs[qn] && p->ZN[KP]<p->B270_ze[qn])
	set_cell(p,a,p->B270_n[qn],p->B270_d50[qn],p->B270_alpha[qn],p->B270_beta[qn]);
    
    // Vertical Cylinder
    for(qn=0;qn<p->B274;++qn)
    LOOP
    {
        double  r = sqrt( pow(p->XP[IP]-p->B274_xc[qn],2.0)+pow(p->YP[JP]-p->B274_yc[qn],2.0));
        
        if(r<=p->B274_r[qn] && p->pos_z()>p->B274_zs[qn] && p->pos_z()<=p->B274_ze[qn])
        set_cell(p,a,p->B274_n[qn],p->B274_d50[qn],p->B274_alpha[qn],p->B274_beta[qn]);
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

		LOOP
		if(p->pos_x()>=xs && p->pos_x()<xe
		&& p->pos_y()>=p->B281_ys[qn] && p->pos_y()<p->B281_ye[qn]
		&& p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_x()-p->B281_xs[qn])+p->B281_zs[qn] )
		set_cell(p,a,p->B281_n[qn],p->B281_d50[qn],p->B281_alpha[qn],p->B281_beta[qn]);
    }
    
    // Wedge y-dir
    for(qn=0;qn<p->B282;++qn)
    {
		zmin=MIN(p->B282_zs[qn],p->B282_ze[qn]);
        
            if(p->B282_ys[qn]<=p->B282_ye[qn])
            {
            ys = p->B282_ys[qn];
            ye = p->B282_ye[qn];
            }
            
            if(p->B282_ys[qn]>p->B282_ye[qn])
            {
            ys = p->B282_ye[qn];
            ye = p->B282_ys[qn];
            }

		slope=(p->B282_ze[qn]-p->B282_zs[qn])/(p->B282_ye[qn]-p->B282_ys[qn]);

		LOOP
		if(p->pos_x()>=p->B282_xs[qn] && p->pos_x()<p->B282_xe[qn]
		&& p->pos_y()>=ys && p->pos_y()<ye
		&& p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_y()-p->B282_ys[qn])+p->B282_zs[qn] )
		set_cell(p,a,p->B282_n[qn],p->B282_d50[qn],p->B282_alpha[qn],p->B282_beta[qn]);
    }
    
    // Plate x-dir
    for(qn=0;qn<p->B291;++qn)
    {
		zmin=MIN(p->B291_zs[qn],p->B291_ze[qn]);
        zmax=MAX(p->B291_zs[qn],p->B291_ze[qn]) + p->B291_d[qn];
        
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

		LOOP
		if(p->pos_x()>=xs && p->pos_x()<xe
		&& p->pos_y()>=p->B291_ys[qn] && p->pos_y()<p->B291_ye[qn]
        
		&& p->pos_z()>=zmin 
        && p->pos_z()<=zmax 
        
        && p->pos_z()<slope*(p->pos_x()-p->B291_xs[qn])+p->B291_zs[qn]+p->B291_d[qn] // upper
        && p->pos_z()>slope*(p->pos_x()-p->B291_xs[qn])+p->B291_zs[qn]) // lower
		set_cell(p,a,p->B291_n[qn],p->B291_d50[qn],p->B291_alpha[qn],p->B291_beta[qn]);
    }
}
