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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

particle_pls::particle_pls(lexer* p, fdm *a, ghostcell* pgc) : norm_vec(p), phimax(p),phimin(p),phiold(p),
							  posnum(p), negnum(p),
                             zero (0.0), epsi(1.5*p->DXM),dx(p->DXM),rmin(0.1*p->DXM),
                             rmax(0.5*p->DXM),pnum(p->F32),ipolval(p->F31), irand(100000), drand(100000.0),
							  nu(1.0e-10*p->DXM)
{
    pcount=0;
    ncount=0;
    posactive=0;
    negactive=0;
	printcount=0;
    
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_PLS",0777);
}

particle_pls::~particle_pls()
{
}

void particle_pls::start(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow)
{ 

	starttime=pgc->timer();
	
	posactive_old=posactive;
	negactive_old=negactive;

    dgc_update(p,a,pgc);
    advect(p,a,pgc,pos,posflag,posactive);
    advect(p,a,pgc,neg,negflag,negactive);
	particlex(p,a,pgc);
    remove(p,a,pgc);
	
	if((p->count%p->F34==0 || p->count==0) && p->F34>0)
	{
    print_vtu(p,a,pgc,pos,posflag,posactive,1);
	print_vtu(p,a,pgc,neg,negflag,negactive,2);
	++printcount;
	}
	
    correct(p,a,pgc,pflow);
	xupdate(p,a,pgc);
	parcount(p,a,pgc);
	random_delete(p,a,pgc);
    reseed(p,a,pgc,0.5);
	setradius(p,a);    
    vel_setback(p,a,pgc);
    pgc->start4(p,a->phi,gcval_phi);

	posbalance = posactive - posactive_old;
	negbalance = negactive - negactive_old;

    gnegactive = pgc->globalisum(negactive);
    gposactive = pgc->globalisum(posactive);
    gpcount = pgc->globalisum(pcount);
    gncount = pgc->globalisum(ncount);
    gcorrected = pgc->globalisum(corrected);
    gremoved = pgc->globalisum(removed);
    greseeded = pgc->globalisum(reseeded);
    gxchange = pgc->globalisum(xchange);
	gposbalance = pgc->globalisum(posbalance);
	gnegbalance = pgc->globalisum(negbalance);
	
	p->plstime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
	{
    cout<<"PLS. neg: "<<gnegactive<<" pos: "<<gposactive<<" n: "<<gncount<<" p: "<<gpcount<<" nbal: "<<gnegbalance<<" pbal: "<<gposbalance<<endl;
	cout<<"CORR: *"<<gcorrected<<"* rem: "<<gremoved<<" res: "<<greseeded<<" X: "<<gxchange<<" | plstime: "<<p->plstime<<endl;
	}
}






