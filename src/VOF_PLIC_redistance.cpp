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

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"


void VOF_PLIC::redistance
(
	fdm* a, 
	lexer* p, 
	convection* pconvec, 
	ghostcell* pgc, 
	ioflow* pflow, 
	int bandWidth
)
{    
    int ip, jp, kp, flag;
    double xp, yp, zp;
    field4 changedFlag(p);
	
	LOOP
	{
		changedFlag(i,j,k) = 0;
	}
	

    //- Main loop
    LOOP
    {
		reconstructPlane(a, p);
		
        if (a->vof(i, j, k) < 0.999 && a->vof(i, j, k) > 0.001)
        {
	//		cout<<"\nCell ID: "<<i<<" "<<k<<" "<<p->pos_x()<<" "<<p->pos_z()<<" "<<a->vof(i, j, k)<<endl;
			
			for (int ii = -bandWidth; ii < (bandWidth + 1); ii++)
			{
				for (int jj = -bandWidth; jj < (bandWidth + 1); jj++)
				{  
					for (int kk = -bandWidth; kk < (bandWidth + 1); kk++)
					{
						ip = i + ii;
						jp = j + jj;
						kp = k + kk;
  
						if 
						(
							ip >= 0 && jp >= 0 && kp >= 0 && 
							ip < p->knox && jp < p->knoy && kp < p->knoz
						)  
						{ 
//cout<<"Cell considered: "<<ip<<" "<<kp<<endl;	

							//- Calculate closest point on cell boundary
							if (i != ip || j != jp || k != kp)
							{
								flag = calcBoundaryPoint(a, p, ip, jp, kp, changedFlag);
							}
							else
							{
								flag = 1;
							}
							
							if (flag == 1)
							{ 
								//- Calculate projection point
								flag = 
									calcProjectionPoint
									(
										a, p, xp, yp, zp, ip, jp, kp, changedFlag
									);
								
                              if (flag == 2)
								{ 
									//- Calculate closest point on boundary of 
									//- interface segment
									//calcSegmentPoint
									(
										a, p, xp, yp, zp, ip, jp, kp, changedFlag
									);
								}                               
                            }
                        }
                    }
                }
            }   
        }
    }
}


void VOF_PLIC::calcSegmentPoint
(
    fdm* a, 
    lexer* p,
	double xp, 
	double yp,
	double zp,
    int ip, 
	int jp, 
	int kp,
	field4& changedFlag
)
{
	double xs, ys, zs;
	
	double xi = p->pos_x();
	double yi = p->pos_y();
	double zi = p->pos_z();
	
	double xoff = max(fabs(xp - xi) - 0.5*p->DXM, 0.0);
	double yoff = max(fabs(yp - yi) - 0.5*p->DXM, 0.0);
	double zoff = max(fabs(zp - zi) - 0.5*p->DXM, 0.0);
                                    
	double xfc = xi + SIGN(xp - xi)*0.5*p->DXM;
	double yfc = yi + SIGN(yp - yi)*0.5*p->DXM;
	double zfc = zi + SIGN(zp - zi)*0.5*p->DXM;
                                    
	zs = zfc;
                                    
	if (xoff*fabs(nx(i, j, k)) >= yoff*fabs(ny(i, j, k)))
	{
		xs = xfc;
		ys = (alpha(i, j, k)*p->DXM - nz(i, j, k)*zfc - nx(i, j, k)*xfc)/ny(i, j, k);
	}
    else 
	{
		ys = yfc;
		xs = (alpha(i, j, k)*p->DXM - nz(i, j, k)*zfc - ny(i, j, k)*yfc)/nx(i, j, k);                                     
    }
                                    
	double dx = xp - xs;
	double dy = yp - ys;
	double dz = zp - zs;
	
	double d = sqrt(dx*dx + dy*dy + dz*dz);
                        
	if (((fabs(d) < fabs(a->phi(ip, jp, kp))) && changedFlag(ip,jp,kp) == 1) || changedFlag(ip,jp,kp) == 0)
	{
		a->phi(ip, jp, kp) = fabs(d)*SIGN(a->vof(ip, jp, kp) - 0.5);    
		changedFlag(ip,jp,kp) = 1;  
	}
}


int VOF_PLIC::calcProjectionPoint
(
    fdm* a, 
    lexer* p, 
	double& xp, 
	double& yp,
	double& zp,
    int ip, 
	int jp, 
	int kp,
	field4& changedFlag
)
{
	double xi = p->pos_x();
	double yi = p->pos_y();
	double zi = p->pos_z();
	
	xp = (double(ip) + 0.5)*p->DXM + p->originx;
	yp = (double(jp) + 0.5)*p->DXM + p->originy;
	zp = (double(kp) + 0.5)*p->DXM + p->originz;
    
	
	double d = 
		calcDistance
		(
			alpha(i, j, k)*p->DXM, 
			xp-(xi - 0.5*p->DXM), 
			yp-(yi - 0.5*p->DXM), 
			zp-(zi - 0.5*p->DXM)
		);
 
	xp += nx(i, j, k)*d;
	yp += ny(i, j, k)*d;
	zp += nz(i, j, k)*d;
                         
	if 
	( 
		xp >= (xi - 0.5*p->DXM) 
        && xp <= (xi + 0.5*p->DXM)
        && yp >= (yi - 0.5*p->DXM) 
        && yp <= (yi + 0.5*p->DXM)
        && zp >= (zi - 0.5*p->DXM) 
        && zp <= (zi + 0.5*p->DXM)
	)
	{
		if (((fabs(d) < fabs(a->phi(ip, jp, kp))) && changedFlag(ip,jp,kp) == 1) || changedFlag(ip,jp,kp) == 0)
		{
			a->phi(ip, jp, kp) = fabs(d)*SIGN(a->vof(ip, jp, kp) - 0.5); 
			changedFlag(ip,jp,kp) = 1; 
		}
     
		return 0;
	}
	else
	{
		return 2;
	}
}


int VOF_PLIC::calcBoundaryPoint
(
    fdm* a, 
    lexer* p, 
    int ip, 
	int jp, 
	int kp,
	field4& changedFlag
)
{
    int l = max(-1, min(1, ip - i));
    int m = max(-1, min(1, jp - j));
    int n = max(-1, min(1, kp - k));
    
	double xv = 0.0;
	if (l == 0) xv = 0.5;
	else if (l == 1) xv = 1.0;
	double yv = 0.0;
	if (m == 0) yv = 0.5;
	else if (m == 1) yv = 1.0;
	double zv = 0.0;
	if (n == 0) zv = 0.5;
	else if (n == 1) zv = 1.0;

    double d =  calcDistance(alpha(i, j, k), xv, yv, zv);
	
	if (d*SIGN(a->vof(ip, jp, kp) - 0.5) > 0.0)
	{
		double dx = 
			fabs
			(
				(double((i + l/2) + 0.5)*p->DXM + p->originx) 
				- (double(ip + 0.5)*p->DXM + p->originx)
			);
		double dy = 
			fabs
			(
				(double((j + m/2) + 0.5)*p->DXM + p->originy) 
				- (double(jp + 0.5)*p->DXM + p->originy)
			);
		double dz = 
			fabs
			(
				(double((k + n/2) + 0.5)*p->DXM + p->originz) 
				- (double(kp + 0.5)*p->DXM + p->originz)
			);
		
		d = sqrt(dx*dx + dy*dy + dz*dz);

		if (((fabs(d) < fabs(a->phi(ip, jp, kp))) && changedFlag(ip,jp,kp) == 1) || changedFlag(ip,jp,kp) == 0)
		{
			a->phi(ip, jp, kp) = fabs(d)*SIGN(a->vof(ip, jp, kp) - 0.5); 
			
			changedFlag(ip,jp,kp) = 1; 
			
			return 0;
		}    
	}
	else
	{
		return 1;
	}
}


double VOF_PLIC::calcDistance(double alpha, double x, double y, double z)
{
    return 
    (
        (
            alpha 
            - (nx(i, j, k)*x + ny(i, j, k)*y + nz(i, j, k)*z)
        )
        /
        (
            sqrt
            (
                  nx(i, j, k)*nx(i, j, k) 
                + ny(i, j, k)*ny(i, j, k) 
                + nz(i, j, k)*nz(i, j, k)
            )
        )
    );
}
