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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

particle_f::particle_f(lexer* p, fdm *a, ghostcell* pgc) : norm_vec(p), active(p),posnum(p), 
                                epsi(1.5*p->DXM), dx(p->DXM), dy(p->DXM), dz(p->DXM),rmin(0.1*p->DXM),
								rmax(0.5*p->DXM), irand(100000), drand(100000.0)
{
    pcount=0;
    posactive=0;
    
    if(p->I40==0)
    {
	printcount=0;
    p->partprinttime=0.0;
    }
    
    // if(p->F50==1)
	// gcval_phi=51;

	// if(p->F50==2)
	// gcval_phi=52;

	// if(p->F50==3)
	// gcval_phi=53;

	// if(p->F50==4)
	// gcval_phi=54;
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_Particle",0777);
}

particle_f::~particle_f()
{
}

void particle_f::start(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow)
{ 

	starttime=pgc->timer();
	
	posactive_old=posactive;

	if (p->count>=p->Q43)
    	advect(p,a,pgc,pos,posflag,posactive);
	particlex(p,a,pgc);
    remove(p,a,pgc);
	
	print_particles(p,a,pgc);
	
    
	xupdate(p,a,pgc);
	parcount(p,a,pgc); 

    pgc->start4(p,a->phi,gcval_phi);

	posbalance = posactive - posactive_old;

    gposactive = pgc->globalisum(posactive);
    gpcount = pgc->globalisum(pcount);
    gcorrected = pgc->globalisum(corrected);
    gremoved = pgc->globalisum(removed);
    greseeded = pgc->globalisum(reseeded);
    gxchange = pgc->globalisum(xchange);
	gposbalance = pgc->globalisum(posbalance);
	
	p->plstime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
	{
    cout<<"PLS. pos: "<<gposactive<<" p: "<<gpcount<<" pbal: "<<gposbalance<<endl;
	cout<<"CORR: *"<<gcorrected<<"* rem: "<<gremoved<<" res: "<<greseeded<<" X: "<<gxchange<<" | plstime: "<<p->plstime<<endl;
	}
}






