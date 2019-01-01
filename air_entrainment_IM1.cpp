/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"air_entrainment_IM1.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_concentration.h"

air_entrainment_IM1::air_entrainment_IM1(lexer* p, fdm* a, ghostcell *pgc) : concentration_io(p,a), air_entrainment_ibc(p), Cn(p)
{
	gcval_concentration=80;
}

air_entrainment_IM1::~air_entrainment_IM1()
{
}

void air_entrainment_IM1::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, turbulence *pturb, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    clearrhs(p,a);
    pconvec->start(p,a,C,4,a->u,a->v,a->w);
	pdiff->idiff_scalar(p,a,pgc,psolv,C,a->visc,1.0,1.0);
	air_entrainment_ibc_start(p,a,pgc,pturb,C);
	timesource(p,a,C);
	
	psolv->start(p,a,pgc,C,a->xvec,a->rhsvec,4,gcval_concentration,p->N43);
	
	pgc->start4(p,C,gcval_concentration);
	p->concentrationtime=pgc->timer()-starttime;
	p->concentrationiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"concentrationiter: "<<p->concentrationiter<<"  concentrationtime: "<<setprecision(3)<<p->concentrationtime<<endl;

}

void air_entrainment_IM1::timesource(lexer* p, fdm* a, field& fn)
{
    int count=0;
    int q;

    LOOP
    {
        a->M.p[count]+= 1.0/PDT;

        a->rhsvec.V[count] += a->L(i,j,k) + Cn(i,j,k)/PDT;

	++count;
    }
}

void air_entrainment_IM1::ttimesave(lexer* p, fdm* a)
{
    LOOP
    Cn(i,j,k)=C(i,j,k);
}

void air_entrainment_IM1::clearrhs(lexer* p, fdm* a)
{

    int count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	a->L(i,j,k)=0.0;
	++count;
    }

}

