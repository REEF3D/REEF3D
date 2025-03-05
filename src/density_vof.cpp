/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"initialize.h"
#include"density_vof.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

density_vof::density_vof(lexer* p) : epsi(p->F45*p->DXM), eps(2.1*p->DXM)
{
    
    
    double psim;
    int count;
    
    if(p->j_dir==0)        
    p->psi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    p->psi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
    
    
    p->psi0=p->psi;
}

density_vof::~density_vof()
{
}

double density_vof::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{  
    if(p->F92==1)   
    {
        double phival, psiro;
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
        psiro = p->psi;
        
        if(phival>psiro)
            H=1.0;

        if(phival<-psiro)
            H=0.0;

        if(fabs(phival)<=psiro)
            H=0.5*(1.0 + phival/(psiro) + (1.0/PI)*sin((PI*phival)/(psiro)));
    
        roval = p->W1*H + p->W3*(1.0-H);
    }
    
    if(p->F92==2)
    {
        if(a->vof(i,j,k)>0.999 && a->vof(i+aa,j+bb,k+cc)>0.999)
            roval=p->W1;
        else if(a->vof(i,j,k)<0.001 && a->vof(i+aa,j+bb,k+cc)<0.001)
            roval=p->W3;
        else
        {
            roval=((a->vof(i,j,k)*p->W1+(1.0-a->vof(i,j,k))*p->W3)
                    + (a->vof(i+aa,j+bb,k+cc)*p->W1+(1.0-a->vof(i+aa,j+bb,k+cc))*p->W3)
                    ) /2.0;
        }
    }

	return roval;
	
}




