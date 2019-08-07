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

#include"density_heat.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"

density_heat::density_heat(lexer* p, heat *& ppheat) : epsi(p->F45*p->dx), eps(2.1*p->dx)
{
    pheat = ppheat;
}

density_heat::~density_heat()
{
}

double density_heat::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi,r,s;
    double heatval;
	
	ii = aa-aa/(fabs(aa)>0?fabs(aa):1);
	jj = bb-bb/(fabs(bb)>0?fabs(bb):1);
	kk = cc-cc/(fabs(cc)>0?fabs(cc):1);
	

	if(p->D32==2)
	{
       
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
        
        heatval = 0.5*(pheat->val(i,j,k) + pheat->val(i+aa,j+bb,k+cc));
        
        if(p->j_dir==0)
        psi = p->F45*(1.0/1.0)*(p->DXN[IP]+p->DZN[KP]);
        
        if(p->j_dir==1)
        psi = p->F45*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
    
        if(phival>psi)
        H=1.0;

        if(phival<-psi)
        H=0.0;

        if(fabs(phival)<=psi)
        H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
        
        roval = (p->W1+heatval*p->C1)*H + (p->W3+heatval*p->C3)*(1.0-H);
	}
	
	// -----
	
	if(p->D32==3)
	roval = 0.5*(a->ro(i+ii,j+jj,k+kk) + a->ro(i+aa,j+bb,k+cc));
	
	// -----
	
	return roval;		
}




