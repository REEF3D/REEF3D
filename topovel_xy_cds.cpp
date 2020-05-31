/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"topo_vel.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bedconc.h"
#include"topo_relax.h"

void topo_vel::topovel_xy_cds(lexer* p,fdm* a, ghostcell *pgc, double& vx, double& vy, double& vz)
{
	double uvel,vvel,uvel1,vvel1,u_abs;
	double sgx,sgy,sgx1,sgy1;
	double dqx,dqy;
	
	vx=0.0;
	vy=0.0;
	vz=0.0;
	 
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{		/*		
        uvel  = a->P(i,j);
        uvel1 = a->P(i-1,j);
        vvel  = a->Q(i,j);
        vvel1 = a->Q(i,j-1);
        
		sgx=fabs(uvel)>1.0e-10?uvel/fabs(uvel):0.0;
		sgy=fabs(vvel)>1.0e-10?vvel/fabs(vvel):0.0;
        sgx1=fabs(uvel1)>1.0e-10?uvel1/fabs(uvel1):0.0;
		sgy1=fabs(vvel1)>1.0e-10?vvel1/fabs(vvel1):0.0;

        dqx = (a->qbx(i,j)*sgx-a->qbx(i-1,j)*sgx1)/(p->DXN[IP]);
        dqy = (a->qby(i,j)*sgy-a->qby(i,j-1)*sgy1)/(p->DYN[JP]);
        */
        uvel  = 0.5*(a->P(i,j)+a->P(i-1,j));
        vvel  = 0.5*(a->Q(i,j)+a->Q(i,j-1));
        
		sgx=fabs(uvel)>1.0e-10?uvel/fabs(uvel):0.0;
		sgy=fabs(vvel)>1.0e-10?vvel/fabs(vvel):0.0;


        dqx = sgx*(a->qbx(i,j)-a->qbx(i-1,j))/(p->DXN[IP]);
        dqy = sgy*(a->qby(i,j)-a->qby(i,j-1))/(p->DYN[JP]);
		
	// Exner equations
    vz =  -prelax->rf(p,a,pgc)*(1.0/(1.0-p->S24))*(dqx + dqy) + ws*(a->conc(i,j,k) - pcb->cbed(p,a,pgc,a->topo)); 
	}
}