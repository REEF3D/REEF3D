/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"ibcsusp.h"
#include"lexer.h"
#include"fdm.h"
#include"turbulence.h"

ibcsusp::ibcsusp(lexer *p, turbulence *pturb) : bedconc(p,pturb), epsi(1.6*p->dx)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    shields=p->S30;
    visc=p->W2;
    kappa=0.4;
    ks=2.5*d50;
    Rstar=(rhosed-rhowat)/rhowat;
}

ibcsusp::~ibcsusp()
{
}

void ibcsusp::ibcsusp_start(lexer* p,fdm* a,ghostcell *pgc,field& conc)
{
	double h,adist,zdist,concval;
	int kmem;
    
    adist = 2.0*d50;

		GC4LOOP
		if(p->gcb4[n][4]==5)
		{
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
		
		kmem=k;
        
        h=a->phi(i,j,k) + p->DZP[KP];
        zdist = p->DZP[KP];

        concval = cbed(p,a,pgc,a->topo)*pow(((h-zdist)/zdist)*(adist/(h-adist)),zdist);
		
		k=kmem;
		conc(i,j,k) =  concval;
		}
}
