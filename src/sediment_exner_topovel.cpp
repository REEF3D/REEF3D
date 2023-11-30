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

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bedconc.h"
#include"topo_relax.h"
#include"sediment_exnerdisc.h"
#include"sediment_fdm.h"

void sediment_exner::topovel(lexer* p, ghostcell *pgc, sediment_fdm *s, double& vx, double& vy, double& vz)
{
	double uvel,vvel,u_abs;
	double signx,signy;
	double dqx,dqy;
    double qx1,qx2,q1x,qy2;
    
    double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
    double sgx1,sgx2,sgy1,sgy2;
    double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
	
	vx=0.0;
	vy=0.0;
	vz=0.0;
	 
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
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
        if(p->S17==0)
        {
        dqx = pdx->sx(p,s->qbe,sgx1,sgx2);
        dqy = pdx->sy(p,s->qbe,sgy1,sgy2);
        }
        
        if(p->S17==1)
        {
        dqx = pdx->sx(p,s->qb,sgx1,sgx2);
        dqy = pdx->sy(p,s->qb,sgy1,sgy2);
        }
        
        vx=dqx;
        vy=dqy;
		
    // Exner equations
        // eq
        if(p->S17==0)
        vz =  -s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(dqx + dqy); //+ susp_qb(p,pgc,s);
        
        
        // non-eq
        if(p->S17==1)
        {
        Ls = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j), 0.0)*d50;
        
        vz =  s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(1.0/(Ls>1.0e-10?Ls:1.0e10))*(s->qb(i,j)-s->qbe(i,j));// + ws*(s->conc(i,j,k) - pcb->cbed(p,pgc,s)); 
        }
	}
    
    /*
    if(p->mpirank==7)
    SLICELOOP4
    if(p->flagslice4[IJ]>0 && p->flagslice4[Ip1J]<0)
    cout<<s->shields_eff(i-2,j)<<" "<<s->shields_eff(i-1,j)<<" "<<s->shields_eff(i,j)<<" "<<s->shields_eff(i+1,j)<<endl;
    */
    
    
}

