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

#include"fixtimestep_ptf.h"
#include<iomanip>
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"turbulence.h"

fixtimestep_ptf::fixtimestep_ptf(lexer* p)
{
}

fixtimestep_ptf::~fixtimestep_ptf()
{
}

void fixtimestep_ptf::start(fdm_ptf *e, lexer* p,ghostcell* pgc, turbulence *pturb)
{
   	p->dt=p->N49;

   	// maximum velocities


	ULOOP
	p->umax=MAX(p->umax,fabs(e->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);


	VLOOP
	p->vmax=MAX(p->vmax,fabs(e->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	WLOOP
	p->wmax=MAX(p->wmax,fabs(e->w(i,j,k)));

	p->wmax=pgc->globalmax(p->wmax);



    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    }

// maximum viscosity
	LOOP
	p->viscmax=MAX(p->viscmax, e->visc(i,j,k)+e->eddyv(i,j,k));

	p->viscmax=pgc->globalmax(p->viscmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"viscmax: "<<p->viscmax<<endl;
	//----kin
	LOOP
	p->kinmax=MAX(p->kinmax,pturb->kinval(i,j,k));

	p->kinmax=pgc->globalmax(p->kinmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"kinmax: "<<p->kinmax<<endl;

	//---eps
    LOOP
	p->epsmax=MAX(p->epsmax,pturb->epsval(i,j,k));

	p->epsmax=pgc->globalmax(p->epsmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"epsmax: "<<p->epsmax<<endl;


	//---press
    LOOP
    {
	p->pressmax=MAX(p->pressmax,e->press(i,j,k));
	p->pressmin=MIN(p->pressmin,e->press(i,j,k));
    }

	p->pressmax=pgc->globalmax(p->pressmax);
	p->pressmin=pgc->globalmin(p->pressmin);

}

void fixtimestep_ptf::ini(fdm_ptf *e, lexer* p,ghostcell* pgc)
{
    p->dt=p->N49;

}
