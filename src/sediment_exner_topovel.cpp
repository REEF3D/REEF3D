/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bedconc_VR.h"
#include"topo_relax.h"
#include"sediment_exnerdisc.h"
#include"sediment_fdm.h"

void sediment_exner::topovel1(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
	double dqx,dqy;
    double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
    double sgx1,sgx2,sgy1,sgy2;
    double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
	

    SLICELOOP4
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{						
        ux1=s->P(i-1,j);
        vx1=0.25*(s->Q(i,j)+s->Q(i-1,j)+s->Q(i,j-1)+s->Q(i-1,j-1)); 
        
        ux2=s->P(i,j);
        vx2=0.25*(s->Q(i,j)+s->Q(i+1,j)+s->Q(i,j-1)+s->Q(i+1,j-1)); 
        
        
        uy1=0.25*(s->P(i,j-1)+s->P(i,j)+s->P(i-1,j-1)+s->P(i-1,j));
        vy1=s->Q(i,j-1); 
        
        uy2=0.25*(s->P(i,j)+s->P(i,j+1)+s->P(i-1,j)+s->P(i-1,j+1));
        vy2=s->Q(i,j); 
        
        
        ux1_abs = sqrt(ux1*ux1 + vx1*vx1);
        ux2_abs = sqrt(ux2*ux2 + vx2*vx2);
        
        uy1_abs = sqrt(uy1*uy1 + vy1*vy1);
        uy2_abs = sqrt(uy2*uy2 + vy2*vy2);
            
        sgx1=fabs(ux1_abs)>1.0e-10?ux1/fabs(ux1_abs):0.0;
        sgx2=fabs(ux2_abs)>1.0e-10?ux2/fabs(ux2_abs):0.0;
        
        sgy1=fabs(uy1_abs)>1.0e-10?vy1/fabs(uy1_abs):0.0;
        sgy2=fabs(uy2_abs)>1.0e-10?vy2/fabs(uy2_abs):0.0;
        

        // complete q
        dqx = pdx->sx(p,s->qb,sgx1,sgx2);
        dqy = pdx->sy(p,s->qb,sgy1,sgy2);
        

    // Exner equations
        s->vz(i,j) =  -s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(dqx + dqy + susp_qb(p,pgc,s));
        
        if(p->DFBED[IJ]<0)
        s->vz(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,s->vz,1);
}

void sediment_exner::topovel2(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
	double dqx,dqy;
    double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
    double sgx1,sgx2,sgy1,sgy2;
    double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
    
    
    double u_abs,signx,signy;
    double uvel,vvel;
    
    SEDSLICELOOP
    {
        pip=1;
        uvel=0.5*(s->P(i,j)+s->P(i-1,j));
        pip=0;

        pip=2;
        vvel=0.5*(s->Q(i,j)+s->Q(i,j-1));
        pip=0;
		
		u_abs = sqrt(uvel*uvel + vvel*vvel);
		signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
		signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;
        
        qbx(i,j) = signx*s->qb(i,j);
        qby(i,j) = signy*s->qb(i,j);
    }
    
    pgc->gcsl_start4(p,qbx,1);
    pgc->gcsl_start4(p,qby,1);
    
	

    SEDSLICELOOP
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{						
        ux1=s->P(i-1,j);
        vx1=0.25*(s->Q(i,j)+s->Q(i-1,j)+s->Q(i,j-1)+s->Q(i-1,j-1)); 
        
        ux2=s->P(i,j);
        vx2=0.25*(s->Q(i,j)+s->Q(i+1,j)+s->Q(i,j-1)+s->Q(i+1,j-1)); 
        
        
        uy1=0.25*(s->P(i,j-1)+s->P(i,j)+s->P(i-1,j-1)+s->P(i-1,j));
        vy1=s->Q(i,j-1); 
        
        uy2=0.25*(s->P(i,j)+s->P(i,j+1)+s->P(i-1,j)+s->P(i-1,j+1));
        vy2=s->Q(i,j); 
        
        
        ux1_abs = sqrt(ux1*ux1 + vx1*vx1);
        ux2_abs = sqrt(ux2*ux2 + vx2*vx2);
        
        uy1_abs = sqrt(uy1*uy1 + vy1*vy1);
        uy2_abs = sqrt(uy2*uy2 + vy2*vy2);
            
        sgx1=fabs(ux1_abs)>1.0e-10?ux1/fabs(ux1_abs):0.0;
        sgx2=fabs(ux2_abs)>1.0e-10?ux2/fabs(ux2_abs):0.0;
        
        sgy1=fabs(uy1_abs)>1.0e-10?vy1/fabs(uy1_abs):0.0;
        sgy2=fabs(uy2_abs)>1.0e-10?vy2/fabs(uy2_abs):0.0;
        
        // complete q
        dqx = pdx->sx(p,qbx,sgx1,sgx2);
        dqy = pdx->sy(p,qby,sgy1,sgy2);
        

    // Exner equations
        s->vz(i,j) =  -s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(dqx + dqy + susp_qb(p,pgc,s));
	}
    
    
    //filter(p,pgc,s->vz,1,3);
    
    SEDSLICELOOP
    vztemp(i,j) = s->vz(i,j);
    
    pgc->gcsl_start4(p,vztemp,1);
    
    SEDSLICELOOP
    s->vz(i,j) = 0.5*vztemp(i,j) + 0.125*(vztemp(i-1,j)+vztemp(i+1,j)+vztemp(i,j-1)+vztemp(i,j+1));
    
    pgc->gcsl_start4(p,s->vz,1);
}


void sediment_exner::filter(lexer *p,ghostcell *pgc, slice &f, int outer_iter, int inner_iter)
{
    slice4 h(p),dh(p); 
	
	for(int qn=0;qn<outer_iter;++qn)
	{
		SLICELOOP4
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
		h(i,j) = f(i,j);
		
		pgc->gcsl_start4(p,h,1);
	
        // predictor
		SLICELOOP4
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
		f(i,j) = 0.5*h(i,j) + 0.125*(h(i-1,j) + h(i+1,j) + h(i,j-1) + h(i,j+1));
		
        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            SLICELOOP4
            if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
            dh(i,j) = h(i,j) - f(i,j);
            
            
            SLICELOOP4
            if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
            dh(i,j) = 0.5*dh(i,j) + 0.125*(dh(i-1,j) + dh(i+1,j) + dh(i,j-1) + dh(i,j+1));
            
            SLICELOOP4
            if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
            f(i,j) += dh(i,j);
		}
    }
}