/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"reduction_parker.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

reduction_parker::reduction_parker(lexer *p) : bedslope(p)
{
	eta = 0.85;
}

reduction_parker::~reduction_parker()
{
}


double reduction_parker::start(lexer *p, fdm * a, ghostcell *pgc)
{
    double r=1.0;
	double r1,r2;

	slope(p,a,pgc,teta,alpha,gamma,phi);

	alpha = fabs(alpha);

	mu = atan(1.0/phi);
	d = (4.0/3.0)*mu*0.85*0.74;
	
	pval = (2.0/(1.0-d))*(d/sqrt(1.0 + tan(alpha)*tan(alpha) + tan(teta)*tan(teta)) + sin(teta)/mu);

	qval = ((1.0+d)/(1.0-d))*(1.0/(1.0 + tan(alpha)*tan(alpha) + tan(teta)*tan(teta)))*(-1.0 + ((tan(alpha)*tan(alpha) + tan(teta)*tan(teta))/mu));

	r1 = -0.5*pval - sqrt(pval*pval*0.25 - qval);
	
	r = -0.5*pval + sqrt(pval*pval*0.25 - qval);
	
	if(r<0.0)
    {
    r = 0.1/(fabs(gamma) + 0.0000001)+0.1;
    
    r = MAX(r,0.0);
    r = MIN(r,1.0);
    }
    
	//r=fabs(1.0/(gamma*(180.0/PI)+ 0.0000001));
    
	
	if(  ((1.0 + tan(alpha)*tan(alpha) + tan(teta)*tan(teta))  < 0.0 || (pval*pval*0.25 - qval) < 0.0))
	{
	//r=fabs(1.0/(gamma*(180.0/PI)+ 0.0000001));
    r = 0.1/(fabs(gamma) + 0.0000001)+0.1;
    
    r = MAX(r,0.0);
    r = MIN(r,1.0);
	//cout<<"RRRRRR: "<<r<<endl;
	}
	//r = 0.1/(fabs(gamma) + 0.0000001)+0.1;

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=1.0;
    
    r = MIN(r,2.0);
	
	//if(p->mpirank>0)
	//cout<<"r: "<<r<<"        teta: "<<teta*(180.0/PI)<<" alpha: "<<alpha*(180.0/PI)<<"   gamma: "<<gamma*(180.0/PI)<<"   phi: "<<phi*(180.0/PI)<<endl;

    return r;
}



