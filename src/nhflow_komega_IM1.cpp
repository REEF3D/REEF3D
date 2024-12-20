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

#include"nhflow_komega_IM1.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_strain.h"
#include"solver.h"
#include"nhflow_diffusion.h"
#include"ioflow.h"
#include"nhflow_scalar_convection.h"

nhflow_komega_IM1::nhflow_komega_IM1(lexer* p, fdm_nhf* d, ghostcell *pgc) : nhflow_komega_func(p,d,pgc)
{
	gcval_kin=20;
	gcval_eps=30;
    
    p->Darray(KN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(EN,p->imax*p->jmax*(p->kmax+2));
}

nhflow_komega_IM1::~nhflow_komega_IM1()
{
}

void nhflow_komega_IM1::start(lexer* p, fdm_nhf* d, ghostcell* pgc, nhflow_scalar_convection* pconvec, nhflow_diffusion* pdiff,solver* psolv, ioflow* pflow, vrans *pvrans)
{
	Pk_update(p,d,pgc);
	wallf_update(p,d,pgc,WALLF);

//kin
    starttime=pgc->timer();
	clearrhs(p,d);
    pconvec->start(p,d,KIN,4,d->U,d->V,d->omegaF);
	pdiff->diff_scalar(p,d,pgc,psolv,KIN,1.0);
	kinsource(p,d,pvrans);
	timesource(p,d,KN);
    bckomega_start(p,d,KIN,EPS,gcval_kin);
    bckin_matrix(p,d,KIN,EPS);
    psolv->startV(p,pgc,KIN,d->rhsvec,d->M,4);
    pgc->start20V(p,KIN,gcval_kin);
    kinupdate(p,d,pgc);
	p->kintime=pgc->timer()-starttime;
	p->kiniter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"kin_iter: "<<p->kiniter<<"  kin_time: "<<setprecision(3)<<p->kintime<<endl;

//omega
    starttime=pgc->timer();
	clearrhs(p,d);
    pconvec->start(p,d,EPS,4,d->U,d->V,d->omegaF);
	pdiff->diff_scalar(p,d,pgc,psolv,EPS,1.0);
	epssource(p,d,pvrans);
	timesource(p,d,EN);
    bcomega_matrix(p,d,KIN,EPS);
	psolv->startV(p,pgc,EPS,d->rhsvec,d->M,4);
	epsfsf(p,d,pgc);
	bckomega_start(p,d,KIN,EPS,gcval_eps);
	pgc->start30V(p,EPS,gcval_eps);
	p->epstime=pgc->timer()-starttime;
	p->epsiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"omega_iter: "<<p->epsiter<<"  omega_time: "<<setprecision(3)<<p->epstime<<endl;

    eddyvisc(p,d,pgc,pvrans);
    pflow->turb_relax_nhflow(p,d,pgc,d->EV);
    pgc->start24V(p,d->EV,24);
}

void nhflow_komega_IM1::ktimesave(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
    LOOP
    KN[IJK]=KIN[IJK]; 
}

void nhflow_komega_IM1::etimesave(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
    LOOP
    EN[IJK]=EPS[IJK]; 
}

void nhflow_komega_IM1::kinupdate(lexer *p, fdm_nhf* d, ghostcell *pgc)
{
    LOOP
    d->KIN[IJK]=KIN[IJK]; 
    
    pgc->start4V(p,d->KIN,gcval_kin);
}

void nhflow_komega_IM1::timesource(lexer* p, fdm_nhf* d, double *FN)
{
    count=0;
    LOOP
    {
        d->M.p[count] += 1.0/DT;

        d->rhsvec.V[count] += d->L[IJK] + FN[IJK]/DT;

	++count;
    }
}

void nhflow_komega_IM1::clearrhs(lexer* p, fdm_nhf *d)
{
    count=0;
    LOOP
    {
    d->rhsvec.V[count]=0.0;
	d->L[IJK]=0.0;

	++count;
    }
}
