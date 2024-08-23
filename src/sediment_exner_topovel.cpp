/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

double sediment_exner::topovel(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
	double dqx,dqy;
    double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
    double sgx1,sgx2,sgy1,sgy2;
    double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
	
	vz=0.0;
	 
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
        dqx = pdx->sx(p,s->qbe,sgx1,sgx2);
        dqy = pdx->sy(p,s->qbe,sgy1,sgy2);

        
    // Exner equations
        vz =  -s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(dqx + dqy + susp_qb(p,pgc,s));
        
	}
    
    return vz;
}

