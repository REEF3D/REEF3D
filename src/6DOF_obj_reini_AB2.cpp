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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"

void sixdof_obj::reini_AB2(lexer* p, fdm* a, ghostcell* pgc, field &f)
{	
    LOOP
	{
    if(p->j_dir==0)
    dt.V[IJK] = p->F43*MIN(p->DXP[IP],p->DZP[KP]);
    
    if(p->j_dir==1)
	dt.V[IJK] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	}
	
	reiniter=10;
	
	
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<endl<<"initializing fb..."<<endl<<endl;
	reiniter=10;
	}

    for(int q=0;q<reiniter;++q)
	{
		prdisc->start(p,a,pgc,f,L,5);

		if(q==0)
		ALOOP
		frk1.V[IJK]=L.V[IJK];


		ALOOP
		{
		f.V[IJK] += dt.V[IJK]*0.5*(3.0*L.V[IJK] - frk1.V[IJK]);

		frk1.V[IJK]=L.V[IJK];
		}

	pgc->start4a(p,f,50);
	}
}





