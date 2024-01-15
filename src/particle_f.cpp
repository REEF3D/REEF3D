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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

particle_f::particle_f(lexer* p, fdm *a, ghostcell* pgc) : norm_vec(p), active_box(p), active_topo(p), irand(100000), drand(irand), PP(10,p->S20,p->S22,p->S24)
{
    
    if(p->I40==0)
    {
		printcount=0;
		p->partprinttime=0.0;
    }
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_Particle",0777);

	if(0==p->mpirank)
	{
		// tracers_obj test1(5,2);
		//particles_obj test2(1,1,1,1,1);
		// cout<<"test1.size:"<<test1.size<<"test2.size:"<<test2.size<<endl;
		// test1.add_obj(test2);
		// cout<<"test1.size"<<test1.size<<endl; // bug with clear
		// delete &test1; // memory problem
		//delete &test2;
	}
}

particle_f::~particle_f()
{
}

void particle_f::start(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow)
{ 
	starttime=pgc->timer();

	if (p->count>=p->Q43)
	{
		if(p->count%100==0)
			posseed_suspended(p,a,pgc);
		advect(p,a,pgc);
		particlex(p,a,pgc);
		remove(p,a,pgc);
	}

	print_particles(p,a,pgc);

	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchange);
	p->plstime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Particles: active: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | particletime: "<<p->plstime<<endl;
}






