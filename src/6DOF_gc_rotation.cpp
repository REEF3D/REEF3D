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
#include<math.h>

void sixdof_gc::rotation(double &p,double &q,double &r,double phi, double theta, double psi)
{
	double a,b,c;
	
	// phi
	a = p-xg;
	
	b = (q-yg)*cos(phi) - (r-zg)*sin(phi); 
	
	c = (q-yg)*sin(phi) + (r-zg)*cos(phi); 
	
	p=a+xg;
	q=b+yg;
	r=c+zg;	
	
	// theta
	a = (p-xg)*cos(theta) + (r-zg)*sin(theta); 
	
	b = q-yg;
	
	c = -(p-xg)*sin(theta) + (r-zg)*cos(theta); 
	
	p=a+xg;
	q=b+yg;
	r=c+zg;	
	
	// psi
	a = (p-xg)*cos(psi) - (q-yg)*sin(psi); 
	
	b = (p-xg)*sin(psi) + (q-yg)*cos(psi);
	
	c = r-zg;
	
	p=a+xg;
	q=b+yg;
	r=c+zg;	
}


void sixdof_gc::rotation_quaternion(double &p,double &q,double &r,double phi, double theta, double psi)
{
	double a,b,c;
	
	// Distance to origin
    a = p-xg;
    b = q-yg;
    c = r-zg;
 
	// Rotation using Goldstein page 603 (but there is wrong result)
    p = a*(cos(psi)*cos(theta)) + b*(cos(theta)*sin(psi)) - c*sin(theta);
    q = a*(cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi)) + b*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) + c*(cos(theta)*sin(phi));
    r = a*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) + b*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)) + c*(cos(phi)*cos(theta));    

	// Moving back
    p += xg;
    q += yg;
    r += zg;	
}


std::vector<double> sixdof_gc::rotation_R
(
    const std::vector<double>& point
)
{
    std::vector<double> trans_point(3,0.0);    
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            trans_point[i] = trans_point[i] + R_[i][j]*point[j];
        }
    }
    
    return trans_point;
}
