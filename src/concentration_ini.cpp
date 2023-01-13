/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITCOUT
ANY WARRANTY; without even the implied warranty of MERCCANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"concentration_io.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fluid_update_fsf_concentration.h"

void concentration_io::ini(lexer* p, fdm *a, ghostcell* pgc,concentration *pconcentration)
{
	pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconcentration);

	
	double dx=p->DXM;
	double r;
	int qn;

    LOOP
	C(i,j,k)=p->C50_2;

	LOOP
	if(p->XN[IP]>=p->C51 && p->XN[IP]<p->C54
	&& p->YN[JP]>=p->C52 && p->YN[JP]<p->C55
	&& p->ZN[KP]>=p->C53 && p->ZN[KP]<p->C56)
	{
	C(i,j,k)=p->C50_1;
	}



	if(p->C57_1>0||p->C57_2>0||p->C57_3>0||p->C57_4>0)
	{
		LOOP
		if(p->C57_1*p->XP[IP]+ p->C57_2*p->YP[JP] + p->C57_3*p->ZP[KP] < p->C57_4)
		C(i,j,k)=p->C50_1;
	}

	if(p->C58_4>0.0)
	{
		LOOP
		{
		r = sqrt( pow(p->XP[IP]-p->C58_1,2.0)+pow(p->YP[JP]-p->C58_2,2.0)+pow(p->ZP[KP]-p->C58_3,2.0));
		if(r<=p->C58_4)
		C(i,j,k)=p->C50_1;
		}
	}
	
	for(qn=0;qn<p->C75;++qn)
    {
		double xp1,zp1,xp2,zp2,xp3,zp3,xp4,zp4,x0,z0;
		double s,ls,alpha;
		double xc,zc;
		double xr,zr;
		double vel;
		
		x0 = p->C75_x[qn];
		z0 = p->C75_z[qn];
		alpha = fabs(p->C75_a[qn]*(PI/180.0));
		s = p->C75_s[qn];
		ls = p->C75_l[qn];
		vel =p->C75_v[qn];
		
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
									C(i,j,k)=p->C50_1;
									
									a->u(i,j,k) = cos(alpha)*vel;
									a->w(i,j,k) = -sin(alpha)*vel;
								}
						}
				}
		}
		}
		
			
			
		}

    pgc->start4(p,C,80);
	pupdate->start(p,a,pgc);

}

double concentration_io::fx(double x1, double z1, double x2, double z2, double x)
{
	double f;
	
	f = ((z2-z1)/(x2-x1))*(x-x1) + z1;
	
	return f;
	
}

double concentration_io::fz(double x1, double z1, double x2, double z2, double z)
{
	double f;
	
	f = ((x2-x1)/(z2-z1))*(z-z1) + x1;
	
	return f;
	
}
