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

#include"multiphase_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reini.h"

void multiphase_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, printer *pprint, convection *pconvec, solver *psolv)
{	
	int istart, iend, jstart, jend, kstart, kend;
    int qn;
	double xc,yc,zc,r;
	
	LOOP
	ls1(i,j,k)=-1.0;
	
	pgc->start4(p,ls1,50);
	
	LOOP
	ls2(i,j,k)=-1.0;
	
	pgc->start4(p,ls2,50);
		
	
// LS1
	if(p->F360>-1.0e20)
	{
    LOOP
    ls1(i,j,k)=p->F360-p->pos_x();
	}	
	
	if(p->F361>-1.0e20)
	{
    LOOP
    ls1(i,j,k)=p->F361-p->pos_y();
	}	
	
	if(p->F362>-1.0e20)
	{
    LOOP
    ls1(i,j,k)=p->F362-p->pos_z();
	}	
	
	
	for(qn=0;qn<p->F369;++qn)
    {
		double xp1,zp1,xp2,zp2,xp3,zp3,xp4,zp4,x0,z0;
		double s,ls,alpha;
		double xc,zc;
		double xr,zr;
		double vel;
		
		x0 = p->F369_x[qn];
		z0 = p->F369_z[qn];
		alpha = fabs(p->F369_a[qn]*(PI/180.0));
		s = p->F369_s[qn];
		ls = p->F369_l[qn];
		vel = p->F369_v[qn];
		
		xp1 = x0;
		zp1 = z0;
		
		xp2 = s * cos(alpha) + x0;
		zp2 = s * sin(alpha) + z0;
		
		xp4 = ls * cos(PI-alpha) + x0;
		zp4 = ls * sin(PI-alpha) + z0;
		
		xp3 = s * cos(alpha) + xp4;
		zp3 = s * sin(alpha) + zp4;
		
		LOOP
		{
		xc = p->pos_x();
		zc = p->pos_z();
		
		// g1 : P1 - P2
		xr = fz(xp1,zp1,xp2,zp2,zc);
		zr = fx(xp1,zp1,xp2,zp2,xc);
		
		if(xc<xr && zc>zr)
		{	
			// g2 : P4 - P3
			xr = fz(xp4,zp4,xp3,zp3,zc);
			zr = fx(xp4,zp4,xp3,zp3,xc);
			
				if(xc>xr && zc<zr)
				{
					// g3 : P3 - P2
					xr = fz(xp3,zp3,xp2,zp2,zc);
					zr = fx(xp3,zp3,xp2,zp2,xc);
					
						if(xc<xr && zc<zr)
						{
							// g4 : P4 - P1
							xr = fz(xp4,zp4,xp1,zp1,zc);
							zr = fx(xp4,zp4,xp1,zp1,xc);
							
								if(xc>xr && zc>zr)
								{
									ls1(i,j,k)=1;
									
									a->u(i,j,k) = cos(alpha)*vel;
									a->w(i,j,k) = -sin(alpha)*vel;
								}
						}
				}
		}
		}
		
			
		}
	
// F370
    for(qn=0;qn<p->F370;++qn)
    {
        istart = p->posc_i(p->F370_xs[qn]);
        iend = p->posc_i(p->F370_xe[qn]);
        
        jstart = p->posc_j(p->F370_ys[qn]);
        jend = p->posc_j(p->F370_ye[qn]);
        
        kstart = p->posc_k(p->F370_zs[qn]);
        kend = p->posc_k(p->F370_ze[qn]);
    
        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        ls1(i,j,k)=1;
    }
	
	for(qn=0;qn<p->F371;++qn)
    {
        istart = p->posc_i(p->F371_xs[qn]);
        iend = p->posc_i(p->F371_xe[qn]);
        
        jstart = p->posc_j(p->F371_ys[qn]);
        jend = p->posc_j(p->F371_ye[qn]);
        
        kstart = p->posc_k(p->F371_zs[qn]);
        kend = p->posc_k(p->F371_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        ls1(i,j,k)=-1.0;
    }
	
// F374	
	for(qn=0;qn<p->F374;++qn)
	{

		xc = p->F374_xc[qn] - p->originx;
		zc = p->F374_zc[qn] - p->originz;

		LOOP
		{
		r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

		if(r<=p->F374_r[qn])
		ls1(i,j,k)=1.0;
		
		if(r>p->F374_r[qn])
		ls1(i,j,k)=-1.0;
		}
	}
	
// F375
	for(qn=0;qn<p->F375;++qn)
	{

		xc = p->F375_xc[qn] - p->originx;
		zc = p->F375_zc[qn] - p->originz;

		LOOP
		{
		r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

		if(r<=p->F375_r[qn])
		ls1(i,j,k)=-1.0;
		
		if(r>p->F375_r[qn])
		ls1(i,j,k)=1.0;
		}
	}
	
// LS2
	if(p->F380>-1.0e20)
	{
    LOOP
    ls2(i,j,k)=p->F380-p->pos_x();
	}	
	
	if(p->F381>-1.0e20)
	{
    LOOP
    ls2(i,j,k)=p->F381-p->pos_y();
	}	
	
	if(p->F382>-1.0e20)
	{
    LOOP
    ls2(i,j,k)=p->F382-p->pos_z();
	}	


    for(qn=0;qn<p->F390;++qn)
    {
        istart = p->posc_i(p->F390_xs[qn]);
        iend = p->posc_i(p->F390_xe[qn]);
        
        jstart = p->posc_j(p->F390_ys[qn]);
        jend = p->posc_j(p->F390_ye[qn]);
        
        kstart = p->posc_k(p->F390_zs[qn]);
        kend = p->posc_k(p->F390_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        ls2(i,j,k)=1;
    }
	
	for(qn=0;qn<p->F391;++qn)
    {
        istart = p->posc_i(p->F391_xs[qn]);
        iend = p->posc_i(p->F391_xe[qn]);
        
        jstart = p->posc_j(p->F391_ys[qn]);
        jend = p->posc_j(p->F391_ye[qn]);
        
        kstart = p->posc_k(p->F391_zs[qn]);
        kend = p->posc_k(p->F391_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        ls2(i,j,k)=-1;
    }
	
	// F394	
	for(qn=0;qn<p->F394;++qn)
	{

		xc = p->F394_xc[qn] - p->originx;
		zc = p->F394_zc[qn] - p->originz;

		LOOP
		{
		r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

		if(r<=p->F394_r[qn])
		ls2(i,j,k)=1.0;
		
		if(r>p->F394_r[qn])
		ls2(i,j,k)=-1.0;
		}
	}
	
// F395
	for(qn=0;qn<p->F395;++qn)
	{

		xc = p->F395_xc[qn] - p->originx;
		zc = p->F395_zc[qn] - p->originz;

		LOOP
		{
		r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

		if(r<=p->F395_r[qn])
		ls2(i,j,k)=-1.0;
		
		if(r>p->F395_r[qn])
		ls2(i,j,k)=1.0;
		}
	}
    
    // F374	
	for(qn=0;qn<p->F374;++qn)
	{

		xc = p->F374_xc[qn] - p->originx;
		zc = p->F374_zc[qn] - p->originz;

		LOOP
		{
		r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

		if(r<=p->F374_r[qn])
		ls2(i,j,k)=-1.0;
		}
	}
	
	pgc->start4(p,ls1,50);
	pgc->start4(p,ls2,50);
	
	preini->start(a,p,ls1, pgc, pflow);
	preini->start(a,p,ls2, pgc, pflow);
	
	update(p,a,pgc);
}

double multiphase_f::fx(double x1, double z1, double x2, double z2, double x)
{
	double f;
	
	f = ((z2-z1)/(x2-x1))*(x-x1) + z1;
	
	return f;
	
}

double multiphase_f::fz(double x1, double z1, double x2, double z2, double z)
{
	double f;
	
	f = ((x2-x1)/(z2-z1))*(z-z1) + x1;
	
	return f;
	
}
