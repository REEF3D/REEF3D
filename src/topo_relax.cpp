/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"topo_relax.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

topo_relax::topo_relax(lexer *p) 
{
	p->Darray(betaS73,p->S73);
	p->Darray(tan_betaS73,p->S73);
	p->Darray(dist_S73,p->S73);


	for(n=0;n<p->S73;++n)
	betaS73[n] = (p->S73_b[n]+90.0)*(PI/180.0);

	for(n=0;n<p->S73;++n)
	tan_betaS73[n] = tan(betaS73[n]);
}

topo_relax::~topo_relax()
{
}

void topo_relax::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    
	double relax,distot,distcount,zhval,qbval,cbval;
	double tauval, shearvelval, shieldsval;
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
			zhval = s->bedzh(i,j);
            qbval = s->qbe(i,j);
            cbval =  s->cbe(i,j);
            tauval = s->tau_eff(i,j);
            shearvelval = s->shearvel_eff(i,j);
            shieldsval = s->shields_eff(i,j);
			s->bedzh(i,j)=0.0;
            s->qbe(i,j)=0.0;
            s->cbe(i,j)=0.0;
            s->tau_eff(i,j)=0.0;
            s->shearvel_eff(i,j)=0.0;
            s->shields_eff(i,j)=0.0;
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
			s->bedzh(i,j) += (1.0-relax)*p->S73_val[n] + relax*zhval;
            s->qbe(i,j) +=  relax*qbval;
            s->cbe(i,j) +=  relax*cbval;
            s->tau_eff(i,j)=relax*tauval;
            s->shearvel_eff(i,j)=relax*shearvelval;
            s->shields_eff(i,j)=relax*shieldsval;
			}
			
			
			if(distcount>1)
			{
			s->bedzh(i,j) += ((1.0-relax)*p->S73_val[n] + relax*zhval) * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
            s->qbe(i,j) +=  relax*qbval * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
            s->cbe(i,j) +=  relax*cbval * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
            s->tau_eff(i,j) +=  relax*tauval * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
            s->shearvel_eff(i,j) +=  relax*shearvelval * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
            s->shields_eff(i,j) +=  relax*shieldsval * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
			}
			
			}
		}
    }
	
}

double topo_relax::rf(lexer *p, ghostcell *pgc)
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
            val=(relax);
                
			if(distcount>1)
            val += (relax) * (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20));
            
            //cout<<p->XP[IP]<<" "<<val<<endl;
			}
		}
        
    return val;
    
}

double topo_relax::r1(lexer *p, double x, double threshold)
{
    double r=0.0;

    x=(threshold-fabs(x))/(fabs(threshold)>1.0e-10?threshold:1.0e20);
    x=MAX(x,0.0);
    

    r = 1.0 - (exp(pow(x,3.5))-1.0)/(exp(1.0)-1.0);

    return r;
}

double topo_relax::distcalc(lexer *p ,double x0, double y0, double tan_beta)
{
	double x1,y1;
	double dist=1.0e20;

	x1 = p->pos_x();
	y1 = p->pos_y();
	
	dist = fabs(y1 - tan_beta*x1 + tan_beta*x0 - y0)/sqrt(pow(tan_beta,2.0)+1.0);
	
	return dist;
}




