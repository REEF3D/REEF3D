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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"reinitopo_AB2.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"

reinitopo_AB2::reinitopo_AB2(lexer* p):gradient(p),f(p),frk1(p),frk2(p),L(p),dt(p)
{
	if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	gcval_initopo=150;


	prdisc = new reinidisc_fsf(p);

    
    time_preproc(p);   
}

reinitopo_AB2::~reinitopo_AB2()
{
}

void reinitopo_AB2::start(lexer* p, fdm* a, ghostcell* pgc,field &f)
{
	reiniter=p->S37;
	gcval=gcval_topo;
	
    pgc->start4a(p,f,gcval);
    
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<endl<<"initializing topo..."<<endl<<endl;
    reiniter=2*int(p->maxlength/(p->F43*p->DXM));
	gcval=gcval_initopo;
	pgc->start4a(p,f,gcval);
    
	ALOOP
	L.V[IJK]=frk1.V[IJK]=0.0;
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

	pgc->start4a(p,f,gcval);
	}		
}

void reinitopo_AB2::step(lexer* p, fdm *a)
{

	reiniter=p->S37;
}

void reinitopo_AB2::time_preproc(lexer* p)
{
	
    n=0;
	ALOOP
	{
	dt.V[IJK]= p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	++n;
	}
}


