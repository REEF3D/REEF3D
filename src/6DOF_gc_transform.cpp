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

void sixdof_gc::transform_vec_ES(double xvec,double yvec,double zvec,double &xout,double &yout,double &zout)
{
    xout = xvec*cos(psi)*cos(theta) + yvec*sin(psi)*cos(theta) - zvec*sin(theta);
	
	yout = xvec*(-sin(psi)*cos(phi) + sin(phi)*sin(theta)*cos(psi)) + yvec*(cos(psi)*cos(phi) + sin(phi)*sin(theta)*sin(psi)) + zvec*sin(phi)*cos(theta);

	zout = xvec*(sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi)) + yvec*(-sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)) + zvec*(cos(theta)*cos(phi));
}

void sixdof_gc::transform_vec_SE(double xvec,double yvec,double zvec,double &xout,double &yout,double &zout)
{
	
	xout = xvec*cos(psi)*cos(theta) + yvec*(-sin(psi)*cos(phi) + cos(psi)*sin(theta)*sin(phi)) + zvec*(sin(psi)*sin(phi) + cos(psi)*cos(phi)*sin(theta));
	
	yout = xvec*sin(psi)*cos(theta) + yvec*(cos(psi)*cos(phi) + sin(phi)*sin(theta)*sin(psi)) + zvec*(-cos(psi)*sin(phi) + sin(theta)*sin(psi)*cos(phi));
	
	zout = xvec*(-sin(theta)) + yvec*cos(theta)*sin(phi) + zvec*cos(theta)*cos(phi);
}

void sixdof_gc::transform_angle_ES(double alpha,double beta,double gamma,double &alpha_out,double &beta_out,double &gamma_out)
{
	alpha_out = alpha - gamma*sin(phi);
	
	beta_out = beta*cos(phi) + gamma*cos(theta)*sin(phi);
	
	gamma_out = -beta*sin(phi) + gamma*cos(theta)*cos(phi);
}


void sixdof_gc::transform_angle_SE(double alpha,double beta,double gamma,double &alpha_out,double &beta_out,double &gamma_out)
{
	double alpha1,beta1,gamma1;
	double phival,phimain;
	double thetaval,thetamain;
	double psival,psimain;
	
	
	phi1 = phi;
	phi2 = 0.0;
	
	theta1 = theta;
	theta2 = 0.0;
	
	psi1 = psi;
	psi2 = 0.0;
		
	// theta check
	thetaval = (theta1)/PI;
	
	thetamain = double(int(thetaval));
	
	thetaval -=thetamain;
	
	if(fabs(thetaval)<0.01 || fabs(thetaval)>0.99)
	{
	theta1 = 0.5*theta;
	theta2 = 0.5*theta;
    
    psi1 = 0.5*psi;
	psi2 = 0.5*psi;
    
    phi1 = 0.5*phi;
	phi2 = 0.5*phi;
	}
	
	alpha1 = alpha + beta*sin(phi1)*(sin(theta1)/cos(theta1)) + gamma*cos(phi1)*(sin(theta1)/cos(theta1));
	
	beta1 = beta*cos(phi1) - gamma*sin(phi1);
	
	gamma1 = beta*(sin(phi1)/cos(theta1)) + gamma*(cos(phi1)/cos(theta1));
	
	
	alpha_out = alpha1 + beta1*sin(phi2)*(sin(theta2)/cos(theta2)) + gamma1*cos(phi2)*(sin(theta2)/cos(theta2));
	
	beta_out = beta1*cos(phi2) - gamma1*sin(phi2);
	
	gamma_out = beta1*(sin(phi2)/cos(theta2)) + gamma1*(cos(phi2)/cos(theta2));

}

