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
#include"slice1.h"
#include"slice2.h"
#include"fnpf_weno.h"
 
sflow_sediment_f::sflow_sediment_f(lexer* p, fdm2D *b) : tau(p),taucr(p),alpha(p),teta(p),gamma(p),phi(p),topovel1(p),topovel2(p),fh(p),red(p)
{
    p->sedtime=0.0;
    
    pdx = new fnpf_weno(p);
    
    
    SLICELOOP4
    {
    topovel2(i,j)=0.0;
    topovel1(i,j)=0.0;
    red(i,j)=1.0;
    }
    
    midphi=p->S81*(PI/180.0);
    delta=p->S82*(PI/180.0);
    
    relax_ini(p,b);
}

sflow_sediment_f::~sflow_sediment_f()
{
}

void sflow_sediment_f::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    
}

void sflow_sediment_f::start(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &topovel)
{
    if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		if(p->S42==1 && p->count%p->S44==0)
		sediment_algorithm(p,b,pgc,P,Q,topovel);
		
		if(p->S42==2 && p->simtime>=p->sedsimtime)
		{
		sediment_algorithm(p,b,pgc,P,Q,topovel);
		p->sedsimtime = p->simtime + p->S46;
		}
		
		if(p->S42==3  && p->simtime/p->wT>=p->sedwavetime)
		{
		sediment_algorithm(p,b,pgc,P,Q,topovel);
		p->sedwavetime = p->simtime/p->wT + p->S48;
		}
	}
}

void sflow_sediment_f::sediment_algorithm(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &topovel)
{
    starttime=pgc->timer();
    
    // bedshear
    bedslope(p,b,pgc,P,Q);
    bedshear(p,b,pgc,P,Q);
    shields(p,b,pgc);

    // bedload
    bedload(p,b,pgc);
	filter(p,b,pgc,b->qb,p->S102,p->S103);
    
    // exner
    exner(p,b,pgc,P,Q,topovel);
	filter(p,b,pgc,b->bed,p->S100,p->S101);
    
    // sandslide
    bedslope(p,b,pgc,P,Q);
    sandslide(p,b,pgc,P,Q);
    
    bedchange_update(p,b,pgc);   
    
    
    if(p->mpirank==0 && p->count>0)
    cout<<"Sediment Timestep: "<<p->dtsed<<"  Total Time: "<<setprecision(7)<<p->sedtime<<endl;

	if(p->mpirank==0)
    cout<<"Sediment Time: "<<setprecision(5)<<pgc->timer()-starttime<<endl;
}

