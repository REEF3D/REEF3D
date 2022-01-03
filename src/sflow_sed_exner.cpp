/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"sediment_fou.h"
#include"sediment_cds.h"
#include"sediment_wenoflux.h"
#include"sediment_weno_hj.h"

void sflow_sediment_f::exner(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &topovel)
{

    double maxvel;
    
    double uvel,vvel,u_abs;
	double signx,signy;
	double dqx,dqy;
    double qx1,qx2,q1x,qy2;
    
    double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
    double sgx1,sgx2,sgy1,sgy2;
    double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
	 
    SLICELOOP4
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{						
        
        pip=1;
        uvel=0.5*(P(i,j)+P(i-1,j));
        pip=0;

        pip=2;
        vvel=0.5*(Q(i,j)+Q(i,j-1));
        pip=0;
		
		u_abs = sqrt(uvel*uvel + vvel*vvel);
		signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
		signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;

    
        ux1=P(i-1,j);
        vx1=0.25*(Q(i,j)+Q(i-1,j)+Q(i,j-1)+Q(i-1,j-1)); 
        
        ux2=P(i,j);
        vx2=0.25*(Q(i,j)+Q(i+1,j)+Q(i,j-1)+Q(i+1,j-1)); 
        
        
        uy1=0.25*(P(i,j-1)+P(i,j)+P(i-1,j-1)+P(i-1,j));
        vy1=Q(i,j-1); 
        
        uy2=0.25*(P(i,j)+P(i,j+1)+P(i-1,j)+P(i-1,j+1));
        vy2=Q(i,j); 
        
        
        ux1_abs = sqrt(ux1*ux1 + vx1*vx1);
        ux2_abs = sqrt(ux2*ux2 + vx2*vx2);
        
        uy1_abs = sqrt(uy1*uy1 + vy1*vy1);
        uy2_abs = sqrt(uy2*uy2 + vy2*vy2);
            
        sgx1=fabs(ux1_abs)>1.0e-10?ux1/fabs(ux1_abs):0.0;
        sgx2=fabs(ux2_abs)>1.0e-10?ux2/fabs(ux2_abs):0.0;
        
        sgy1=fabs(uy1_abs)>1.0e-10?vy1/fabs(uy1_abs):0.0;
        sgy2=fabs(uy2_abs)>1.0e-10?vy2/fabs(uy2_abs):0.0;
	


    dqx = pdx->sx(p,b->qb,sgx1,sgx2);
    dqy = pdx->sy(p,b->qb,sgy1,sgy2);
		
	// Exner equation
    topovel(i,j) =  -rf(p,b,pgc)*(1.0/(1.0-p->S24))*(dqx + dqy); 
	}
    
    
    // timestep
	SLICELOOP4
	p->maxtopovel = MAX(fabs(topovel(i,j)),p->maxtopovel);	
    
    p->maxtopovel = pgc->globalmax(p->maxtopovel);
    
        // timestep calculation
        if(p->S15==0)
        p->dtsed=MIN(p->S13, (p->S14*p->DXM)/(fabs(p->maxtopovel)>1.0e-15?p->maxtopovel:1.0e-15));

        if(p->S15==1)
        p->dtsed=MIN(p->dt, (p->S14*p->DXM)/(fabs(p->maxtopovel)>1.0e-15?p->maxtopovel:1.0e-15));
        
        if(p->S15==2)
        p->dtsed=p->S13;

        
        p->dtsed=pgc->timesync(p->dtsed);
        
        p->sedtime+=p->dtsed;
        
        p->maxtopovel=0.0;
        
        // bedchange

        SLICELOOP4
        b->bed(i,j) += p->dtsed*topovel(i,j);


    pgc->gcsl_start4(p,b->bed,50);
}

/*
void sflow_sediment_f::exner(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &topovel)
{
    double uvel,vvel,u_abs;
	double signx,signy;
	double dqx,dqy;
    double maxvel;
	 
    SLICELOOP4
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{						
        pip=1;
        uvel=0.5*(P(i,j)+P(i-1,j));
        pip=0;

        pip=2;
        vvel=0.5*(Q(i,j)+Q(i,j-1));
        pip=0;
		
		u_abs = sqrt(uvel*uvel + vvel*vvel);
		signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
		signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;
	

    dqx = pdx->sx(p,b->qb,signx);
    dqy = pdx->sy(p,b->qb,signy);
		
	// Exner equation
    topovel(i,j) =  -rf(p,b,pgc)*(1.0/(1.0-p->S24))*(dqx*signx + dqy*signy); 
    //b->test(i,j) = topovel(i,j);
	}
    
    
    // timestep
	SLICELOOP4
	p->maxtopovel = MAX(fabs(topovel(i,j)),p->maxtopovel);	
    
    p->maxtopovel = pgc->globalmax(p->maxtopovel);
    
        // timestep calculation
        if(p->S15==0)
        p->dtsed=MIN(p->S13, (p->S14*p->DXM)/(fabs(p->maxtopovel)>1.0e-15?p->maxtopovel:1.0e-15));

        if(p->S15==1)
        p->dtsed=MIN(p->dt, (p->S14*p->DXM)/(fabs(p->maxtopovel)>1.0e-15?p->maxtopovel:1.0e-15));
        
        if(p->S15==2)
        p->dtsed=p->S13;

        
        p->dtsed=pgc->timesync(p->dtsed);
        
        p->sedtime+=p->dtsed;
        
        p->maxtopovel=0.0;
        
        // bedchange

        SLICELOOP4
        b->bed(i,j) += p->dtsed*topovel(i,j);


    pgc->gcsl_start4(p,b->bed,50);
}
*/