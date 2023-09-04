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

#include"density_sf.h"
#include"lexer.h"
#include"fdm.h"

density_sf::density_sf(lexer* p) : epsi(p->F45*p->DXM), eps(2.1*p->DXM)
{
}

density_sf::~density_sf()
{
}

double density_sf::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
    
    if(a->solid(i+aa,j+bb,k+cc)>=0.0 && a->topo(i+aa,j+bb,k+cc)>=0.0)
    phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
    
    if(a->solid(i,j,k)>=0.0 && a->topo(i,j,k)>=0.0)
    if(a->solid(i+aa,j+bb,k+cc)<0.0 || a->topo(i+aa,j+bb,k+cc)<0.0)
    phival = a->phi(i,j,k);
    
    if(a->solid(i+aa,j+bb,k+cc)>=0.0 && a->topo(i+aa,j+bb,k+cc)>=0.0)
    if(a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
    phival = a->phi(i+aa,j+bb,k+cc);

    if(phival>p->psi)
    H=1.0;

    if(phival<-p->psi)
    H=0.0;

    if(fabs(phival)<=p->psi)
    H=0.5*(1.0 + phival/(p->psi) + (1.0/PI)*sin((PI*phival)/(p->psi)));
    


    roval = p->W1*H + p->W3*(1.0-H);

	return roval;		
}




