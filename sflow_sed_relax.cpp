/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

void sflow_sediment_f::relax_ini(lexer *p, fdm2D *b) 
{
	p->Darray(betaS73,p->S73);
	p->Darray(tan_betaS73,p->S73);
	p->Darray(dist_S73,p->S73);


	for(n=0;n<p->S73;++n)
	betaS73[n] = (p->S73_b[n]+90.0)*(PI/180.0);

	for(n=0;n<p->S73;++n)
	tan_betaS73[n] = tan(betaS73[n]);
}

void sflow_sediment_f::relax(lexer *p, fdm2D *b, ghostcell *pgc)
{
	double relax,distot,distcount,zhval;
	
	if(p->S73>0)
	SLICELOOP4
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    {
		distot = 0.0;
		distcount=0;
		for(n=0;n<p->S73;++n)
		{
		dist_S73[n] =  distcalc(p,p->S73_x[n],p->S73_y[n],tan_betaS73[n]);
		
			if(dist_S73[n]<p->S73_dist[n])
			{
			zhval = b->bed(i,j);
			b->bed(i,j)=0.0;
			distot += dist_S73[n];
			++distcount;
			}
		}
		
		for(n=0;n<p->S73;++n)
		{
            if(dist_S73[n]<p->S73_dist[n])
			{
			relax = r1(p,dist_S73[n],p->S73_dist[n]);
			
			if(distcount==1)
			{
			b->bed(i,j) += (1.0-relax)*p->S73_val[n] + relax*zhval;
			}
			
			
			if(distcount>1)
			{
			b->bed(i,j) += ((1.0-relax)*p->S73_val[n] + relax*zhval) * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
			}
			
			}
		}
    }
}

double sflow_sediment_f::rf(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double relax,distot,distcount;
    double val=1.0;
    
        distot = 0.0;
		distcount=0;
		for(n=0;n<p->S73;++n)
		{
		dist_S73[n] =  distcalc(p,p->S73_x[n],p->S73_y[n],tan_betaS73[n]);
		
			if(dist_S73[n]<p->S73_dist[n])
			{
			distot += dist_S73[n];
			++distcount;
			}
		}
		
		for(n=0;n<p->S73;++n)
		{
            if(dist_S73[n]<p->S73_dist[n])
			{
            val=0.0;
			relax = r1(p,dist_S73[n],p->S73_dist[n]);
			
			if(distcount==1)
            val=(1.0-relax);
                
			if(distcount>1)
            val += (1.0-relax) * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
			}
		}
    
    return val;
}

double sflow_sediment_f::r1(lexer *p, double x, double threshold)
{
    double r=0.0;

    x=1.0-x/(fabs(threshold)>1.0e-10?threshold:1.0e20);
    x=MAX(x,0.0);
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(exp(1.0)-1.0);

    return r;
}

double sflow_sediment_f::distcalc(lexer *p, double x0, double y0, double tan_beta)
{
	double x1,y1;
	double dist=1.0e20;

	x1 = p->pos_x();
	y1 = p->pos_y();
	
	dist = fabs(y1 - tan_beta*x1 + tan_beta*x0 - y0)/sqrt(pow(tan_beta,2.0)+1.0);
	
	return dist;
}


