/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"partres2.h"
#include"particles_obj.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

void partres2::relax_ini(lexer *p) 
{
	p->Darray(betaQ73,p->Q73);
	p->Darray(tan_betaQ73,p->Q73);
	p->Darray(dist_Q73,p->Q73);


	for(n=0;n<p->Q73;++n)
	betaQ73[n] = (p->Q73_b[n]+90.0)*(PI/180.0);

	for(n=0;n<p->Q73;++n)
	tan_betaQ73[n] = tan(betaQ73[n]);
}

void partres2::relax(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
}

double partres2::rf(lexer *p, double x1, double y1)
{
    double relax,distot,distcount;
    double val=1.0;
    
        distot = 0.0;
		distcount=0;
		for(n=0;n<p->Q73;++n)
		{
		dist_Q73[n] =  distcalc(p,p->Q73_x[n],p->Q73_y[n],x1,y1,tan_betaQ73[n]);
		
			if(dist_Q73[n]<p->Q73_dist[n])
			{
			distot += dist_Q73[n];
			++distcount;
			}
		}
		
		
		for(n=0;n<p->Q73;++n)
		{
            if(dist_Q73[n]<p->Q73_dist[n])
			{
            val=0.0;
			relax = r1(p,dist_Q73[n],p->Q73_dist[n]);
			
			if(distcount==1)
            val=(relax);
                
			if(distcount>1)
            val += (relax) * (1.0 - dist_Q73[n]/(distot>1.0e-10?distot:1.0e20));
            
            //cout<<p->XP[IP]<<" "<<val<<endl;
			}
		}
        
    return val;
    
}

double partres2::r1(lexer *p, double x, double threshold)
{
    double r=0.0;

    x=(threshold-fabs(x))/(fabs(threshold)>1.0e-10?threshold:1.0e20);
    x=MAX(x,0.0);
    

    r = 1.0 - (exp(pow(x,3.5))-1.0)/(exp(1.0)-1.0);

    return r;
}

double partres2::distcalc(lexer *p ,double x0, double y0 ,double x1, double y1, double tan_beta)
{
	double dist=1.0e20;

	dist = fabs(y1 - tan_beta*x1 + tan_beta*x0 - y0)/sqrt(pow(tan_beta,2.0)+1.0);
	
	return dist;
}




