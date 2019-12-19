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

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"

void sflow_sediment_f::shields(lexer *p, fdm2D *b, ghostcell *pgc)
{
 
    if(p->S80==1)
    parker(p,b,pgc);
    
    if(p->S80==2)
    dey_emp(p,b,pgc);
    
    if(p->S80==3)
    dey_ana(p,b,pgc);
    
    if(p->S80==4)
    fredsoe_long(p,b,pgc);
    
    SLICELOOP4
    {
    taucr(i,j) = (p->S30*fabs(p->W22)*(p->S22-p->W1))*p->S20*red(i,j);
    //b->test(i,j) = red(i,j);
    }
    
    pgc->gcsl_start4(p,b->test,1);
}

void sflow_sediment_f::parker(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double r=1.0;
    double eta = 0.85;
    double mu,d,qval,pval,r1,r2;

    SLICELOOP4
    {
	alpha(i,j) = fabs(alpha(i,j));

	mu = atan(1.0/phi(i,j));
	d = (4.0/3.0)*mu*0.85*0.74;
	
	pval = (2.0/(1.0-d))*(d/sqrt(1.0 + tan(alpha(i,j))*tan(alpha(i,j)) + tan(teta(i,j))*tan(teta(i,j))) + sin(teta(i,j))/mu);

	qval = ((1.0+d)/(1.0-d))*(1.0/(1.0 + tan(alpha(i,j))*tan(alpha(i,j)) + tan(teta(i,j))*tan(teta(i,j))))*(-1.0 + ((tan(alpha(i,j))*tan(alpha(i,j)) + tan(teta(i,j))*tan(teta(i,j)))/mu));

	r1 = -0.5*pval - sqrt(pval*pval*0.25 - qval);
	
	r = -0.5*pval + sqrt(pval*pval*0.25 - qval);
	
	if(r<0.0)
    {
    //r = 0.1/(fabs(gamma(i,j)) + 0.0000001)+0.1;
    
    r = cos(teta(i,j))*(1.0 - tan(teta(i,j)/tan(phi(i,j))));
    r*= cos(alpha(i,j))*(1.0 - pow(tan(alpha(i,j)),2.0)/pow(tan(phi(i,j)),2.0));
    
    r = MAX(r,0.1);
    r = MIN(r,1.0);
    }

	
	if(  ((1.0 + tan(alpha(i,j))*tan(alpha(i,j)) + tan(teta(i,j))*tan(teta(i,j)))  < 0.0 || (pval*pval*0.25 - qval) < 0.0))
	{
    r = 0.1/(fabs(gamma(i,j)) + 0.0000001)+0.1;
    
    
    /*r = cos(teta(i,j))*(1.0 - tan(teta(i,j))/tan(phi(i,j)));
    r*= cos(alpha(i,j))*(1.0 - pow(tan(alpha(i,j)),2.0)/pow(tan(phi(i,j)),2.0));
    
    r = MAX(r,0.1);
    r = MIN(r,2.0);*/
	}

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=1.0;
    

    r=MIN(r,2.0);
    r=MAX(r,0.01);
	
    red(i,j) = r;
    //cout<<r<<endl;
    }
}

void sflow_sediment_f::dey_emp(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double r=1.0;
    double eta = 0.85;

    SLICELOOP4
    {
	alpha(i,j) = fabs(alpha(i,j));

	r = 0.954*pow(1.0-teta(i,j)/phi(i,j), 0.745)*pow(1.0-alpha(i,j)/phi(i,j),0.372);

	if( 1.0-teta(i,j)/phi(i,j) < 0.0 || 1.0-alpha(i,j)/phi(i,j)< 0.0)
    {
    //r = cos(teta(i,j))*(1.0 - tan(teta(i,j))/tan(phi(i,j)));
    //r*= cos(alpha(i,j))*(1.0 - pow(tan(alpha(i,j)),2.0)/pow(tan(phi(i,j)),2.0));
    
	r = 0.954;//0.5/(fabs(gamma(i,j)) + 0.0000001)+0.1;
    }


    if(r<0.0)
    r = 0.01;

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
    
    r=MIN(r,2.0);
    r=MAX(r,0.01);
    
    red(i,j) = r;
    //cout<<r<<" "<<alpha(i,j)<<" "<<teta(i,j)<<endl;
    }
}

void sflow_sediment_f::dey_ana(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double r=1.0;
    double eta = 0.85;

    SLICELOOP4
    {
	alpha(i,j) = fabs(alpha(i,j));

	r = (1.0/((1-eta*tan(phi(i,j)))*tan(phi(i,j))))*( -sin(teta(i,j))  - eta*tan(phi(i,j))*tan(phi(i,j)) * sqrt(cos(teta(i,j))*cos(teta(i,j))-sin(alpha(i,j))*sin(alpha(i,j)))
		+ pow((pow(( sin(teta(i,j)) + eta*tan(phi(i,j))*tan(phi(i,j))*sqrt(cos(teta(i,j))*cos(teta(i,j))-sin(alpha(i,j))*sin(alpha(i,j)))),2.0) +(1.0 - eta*eta*tan(phi(i,j))*tan(phi(i,j)))
		*(cos(teta(i,j))*cos(teta(i,j))*tan(phi(i,j))*tan(phi(i,j)) - sin(alpha(i,j))*sin(alpha(i,j))*tan(phi(i,j))*tan(phi(i,j)) - sin(teta(i,j))*sin(teta(i,j)) - sin(alpha(i,j))*sin(alpha(i,j)) ) ),0.5 ));



	if(  (pow(( sin(teta(i,j)) + eta*tan(phi(i,j))*tan(phi(i,j))*sqrt(cos(teta(i,j))*cos(teta(i,j))-sin(alpha(i,j))*sin(alpha(i,j)))),2.0) +(1.0 - eta*eta*tan(phi(i,j))*tan(phi(i,j)))
		*(cos(teta(i,j))*cos(teta(i,j))*tan(phi(i,j))*tan(phi(i,j)) - sin(alpha(i,j))*sin(alpha(i,j))*tan(phi(i,j))*tan(phi(i,j)) - sin(teta(i,j))*sin(teta(i,j)) - sin(alpha(i,j))*sin(alpha(i,j)) ) )  < 0.0 || cos(teta(i,j))*cos(teta(i,j))-sin(alpha(i,j))*sin(alpha(i,j)) < 0.0)
    {
        r = 0.1/(fabs(gamma(i,j)) + 0.0000001)+0.1;
        
        
        /*r = cos(teta(i,j))*(1.0 - tan(teta(i,j))/tan(phi(i,j)));
        r*= cos(alpha(i,j))*(1.0 - pow(tan(alpha(i,j)),2.0)/pow(tan(phi(i,j)),2.0));
        
        r = MAX(r,0.01);
        r = MIN(r,2.0);*/
    }


    if(r<=0.0)
    r = 0.0001;
    
    r = MIN(r,2.0);

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
    
    r=MIN(r,2.0);
    r=MAX(r,0.01);
    
    red(i,j) = r;
    }
}

void sflow_sediment_f::fredsoe_long(lexer *p, fdm2D *b, ghostcell *pgc)
{  
    double r=1.0;
    
    SLICELOOP4
    {
    r = cos(teta(i,j))*(1.0 - tan(teta(i,j))/tan(phi(i,j)));
    
    r*= cos(alpha(i,j))*(1.0 - pow(tan(alpha(i,j)),2.0)/pow(tan(phi(i,j)),2.0));
    
    r=MIN(r,1.25);
    r=MAX(r,0.01);
    
    red(i,j)=r;
    }
    
}


