/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_bcmom.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"turbulence.h"

nhflow_bcmom::nhflow_bcmom(lexer* p):roughness(p),kappa(0.4)
{
}

nhflow_bcmom::~nhflow_bcmom()
{
}

void nhflow_bcmom::nhflow_bcmom_start(fdm* a, lexer* p,ghostcell *pgc, turbulence *pturb,field& b,int gcval)
{
}

void nhflow_bcmom::roughness_u(lexer* p, fdm_nhf *d, double *U, double *F, slice &WL)
{
    if(p->A519==1)
    {
	k=0;
    
    SLICELOOP4
    {
    dist=p->DZN[KP]*WL(i,j);
	
	ks=d->ks(i,j);


		if(30.0*dist<ks)
		dist=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/ks));

	
	F[IJK] -= (fabs(U[IJK])*U[IJK]*WL(i,j))/(uplus*uplus*dist);
    }
    }
}

void nhflow_bcmom::roughness_v(lexer* p, fdm_nhf *d, double *V, double *G, slice &WL)
{
    if(p->A519==1)
    {
    k=0;
    
    SLICELOOP4
    {
    dist=p->DZN[KP]*WL(i,j);
	
	ks=d->ks(i,j);


		if(30.0*dist<ks)
		dist=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/ks));

	
	G[IJK] -= (fabs(V[IJK])*V[IJK]*WL(i,j))/(uplus*uplus*dist);
    }
    }
	
}

void nhflow_bcmom::roughness_w(lexer* p, fdm_nhf *d, double *W, double *H, slice &WL)
{

}




