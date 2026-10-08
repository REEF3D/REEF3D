/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
    
    p->Darray(dist_S75,p->S75);


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
    // S 73 zones: bed relaxed towards the zone level S73_val, transport quantities damped with r.
    // Overlapping zones: blend weights w_n = (1 - d_n/sum(d))/(N-1), which sum to 1 for any N
    // (the unnormalised weights summed to N-1), and the cell values are taken once before the
    // zone loop (they used to be re-read after the first zone had set them to zero).
	double relax,distot,wn,wsum;
	double zhval,qbval,cbval,tauval,shearvelval,shieldsval;
    double zhnew,qbnew,cbnew,taunew,shearvelnew,shieldsnew;
    int distcount;
    
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
			distot += dist_S73[n];
			++distcount;
			}
		}
        
        if(distcount>0)
        {
        zhval = s->bedzh(i,j);
        qbval = s->qbe(i,j);
        cbval = s->cbe(i,j);
        tauval = s->tau_eff(i,j);
        shearvelval = s->shearvel_eff(i,j);
        shieldsval = s->shields_eff(i,j);
        
        zhnew=qbnew=cbnew=taunew=shearvelnew=shieldsnew=0.0;
        wsum=0.0;
		
            for(n=0;n<p->S73;++n)
            if(dist_S73[n]<p->S73_dist[n])
            {
            relax = r1(p,dist_S73[n],p->S73_dist[n]);
            
            wn = 1.0;
            
            if(distcount>1)
            wn = (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20))/double(distcount-1);
            
            // all zones on top of each other (distot = 0): equal weights
            if(distcount>1 && distot<=1.0e-10)
            wn = 1.0/double(distcount);
            
            wsum += wn;
            
            zhnew += wn*((1.0-relax)*p->S73_val[n] + relax*zhval);
            qbnew += wn*relax*qbval;
            cbnew += wn*relax*cbval;
            taunew += wn*relax*tauval;
            shearvelnew += wn*relax*shearvelval;
            shieldsnew += wn*relax*shieldsval;
            }
        
        wsum = wsum>1.0e-20?wsum:1.0;
        
        s->bedzh(i,j) = zhnew/wsum;
        s->qbe(i,j) = qbnew/wsum;
        s->cbe(i,j) = cbnew/wsum;
        s->tau_eff(i,j) = taunew/wsum;
        s->shearvel_eff(i,j) = shearvelnew/wsum;
        s->shields_eff(i,j) = shieldsnew/wsum;
        }
    }
    
    
    if(p->S75>0)
	SLICELOOP4
    {
		for(n=0;n<p->S75;++n)
		dist_S75[n] =  fabs(p->S75_x[n] - p->pos_x());
		
		for(n=0;n<p->S75;++n)
		{
            if(dist_S75[n]<p->S75_dist[n])
            {
            relax = r1(p,dist_S75[n],p->S75_dist[n]);

            s->bedzh(i,j) = (1.0-relax)*s->bedzh0(i,j) + relax*s->bedzh(i,j);
            s->qbe(i,j) *=  relax;
            s->cbe(i,j) *=  relax;
            s->cb(i,j) *=  relax;
            s->tau_eff(i,j) *=  relax;
            s->shearvel_eff(i,j) *=  relax;
            s->shields_eff(i,j) *=  relax;
            }
		}
    }
	
}

double topo_relax::rf(lexer *p, ghostcell *pgc)
{
    // Exner rate factor, same zone weights as start() (val was reset inside the zone loop,
    // so only the last overlapping zone counted)
    double relax,distot,wn,wsum;
    double val=1.0;
    int distcount;
    
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
		
        if(distcount>0)
        {
        val=0.0;
        wsum=0.0;
        
            for(n=0;n<p->S73;++n)
            if(dist_S73[n]<p->S73_dist[n])
            {
            relax = r1(p,dist_S73[n],p->S73_dist[n]);
            
            wn = 1.0;
            
            if(distcount>1)
            wn = (1.0 - dist_S73[n]/(distot>1.0e-10?distot:1.0e20))/double(distcount-1);
            
            if(distcount>1 && distot<=1.0e-10)
            wn = 1.0/double(distcount);
            
            wsum += wn;
            val += wn*relax;
            }
        
        val = val/(wsum>1.0e-20?wsum:1.0);
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




