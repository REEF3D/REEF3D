/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"simple.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"ioflow.h"

simple::simple(lexer* p,fdm *a) : density(p),apu(p),apv(p),apw(p),pcorr(p)
{
    ini_bcval(p,a);
}

simple::~simple()
{
}

void simple::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{	
    pgc->start1(p,apu,14);
    pgc->start2(p,apv,15);
    pgc->start3(p,apw,16);
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);

    rhs(p,a,pgc,a->u,a->v,a->w,alpha);

    ppois->istart(p,a,apu,apv,apw,pcorr);

	LOOP
    pcorr(i,j,k)=0.0;
    pgc->start4(p,pcorr,gcval_pcorr);

        starttime=pgc->timer();
	psolv->start(p,a,pgc,pcorr,a->xvec,a->rhsvec,5,gcval_press,p->N44);
        endtime=pgc->timer();
    presscorr(p,a,pcorr);
	pgc->start4(p,a->press,gcval_press);
	pgc->start4(p,pcorr,gcval_pcorr);

	ucorr(a,p,pcorr);
	vcorr(a,p,pcorr);
	wcorr(a,p,pcorr);

    convergence(p,a,pgc);

    p->poissoniter+=p->solveriter;
	p->poissontime+=endtime-starttime;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	{
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
	p->poissoniter=0;
	p->poissontime=0.0;
	}

	if(p->mpirank==0 && innercounter==p->N50-1)
	cout<<"pc: "<<p->pcnorm<<" uc: "<<p->ucnorm<<" vc: "<<p->vcnorm<<" wc: "<<p->wcnorm<<endl;
}

void simple::ucorr(fdm* a, lexer* p, field& pcorr)
{
    p->ucnorm=0.0;
	
	ULOOP
	{
	a->u(i,j,k) -=uc=(pcorr(i+1,j,k)-pcorr(i,j,k))
	/(p->DXP[IP]*apu(i,j,k)*roface(p,a,1,0,0));
	p->ucnorm+=uc*uc;
	}
}

void simple::vcorr(fdm* a, lexer* p, field& pcorr)
{
    p->vcnorm=0.0;
	
	VLOOP
	{
	a->v(i,j,k) -=vc= (pcorr(i,j+1,k)-pcorr(i,j,k))
	/(p->DYP[JP]*apv(i,j,k)*roface(p,a,0,1,0));
	p->vcnorm+=vc*vc;
	}
}

void simple::wcorr(fdm* a, lexer* p, field& pcorr)
{
    p->wcnorm=0.0;

	WLOOP
	{
	a->w(i,j,k) -=wc= (pcorr(i,j,k+1)-pcorr(i,j,k))
	/(p->DZP[KP]*apw(i,j,k)*roface(p,a,0,0,1));
	p->wcnorm+=wc*wc;
	}
}

void simple::presscorr(lexer* p, fdm* a, field& pcorr)
{
    p->pcnorm=0.0;
	
    LOOP
    {
    a->press(i,j,k)+=pc=p->N53*pcorr(i,j,k);
    p->pcnorm+=pc*pc;
    }
}

void simple::rhs(lexer *p, fdm* a, ghostcell *pgc, field& u, field& v, field& w, double alpha)
{
    pip=p->Y50;
    
    count=0;
    LOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(p->DXN[IP])
						   -(v(i,j,k)-v(i,j-1,k))/(p->DYN[JP])
						   -(w(i,j,k)-w(i,j,k-1))/(p->DZN[KP]);
                           
    ++count;
    }
    
    pip=0;
}

void simple::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}

void simple::upgrad(lexer*p,fdm* a)
{	
    ULOOP
	a->F(i,j,k)-=PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*roface(p,a,1,0,0));
}

void simple::vpgrad(lexer*p,fdm* a)
{
    VLOOP
	a->G(i,j,k)-=PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))/(p->DYP[JP]*roface(p,a,0,1,0));
}

void simple::wpgrad(lexer*p,fdm* a)
{
    WLOOP
	a->H(i,j,k)-=PORVAL3*(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*roface(p,a,0,0,1));
}

void simple::fillapu(lexer*p,fdm* a)
{
    count=0;
    ULOOP
    {
        apu(i,j,k)=a->M.p[count];

	++count;
    }
}

void simple::fillapv(lexer*p,fdm* a)
{
    count=0;
    VLOOP
    {
        apv(i,j,k)=a->M.p[count];

	++count;
    }
}

void simple::fillapw(lexer*p,fdm* a)
{
    count=0;
    WLOOP
    {
        apw(i,j,k)=a->M.p[count];

	++count;
    }
}

void simple::convergence(lexer *p, fdm *a, ghostcell *pgc)
{
    p->pcnorm = sqrt(pgc->globalsum(p->pcnorm)/double(p->cellnumtot))/(fabs(p->pressmax-p->pressmin)>0.0?fabs(p->pressmax-p->pressmin):1.0);

    p->ucnorm = sqrt(pgc->globalsum(p->ucnorm)/double(p->cellnumtot));
	p->vcnorm = sqrt(pgc->globalsum(p->vcnorm)/double(p->cellnumtot));
	p->wcnorm = sqrt(pgc->globalsum(p->wcnorm)/double(p->cellnumtot));
}

void simple::ini_bcval(lexer *p, fdm *a)
{
	if(p->B76==0)
    gcval_press=40;

    if(p->B76==1)
    gcval_press=41;

    if(p->B76==2)
    gcval_press=42;

    if(p->B76==3)
    gcval_press=43;
	
	if(p->B76==4)
    gcval_press=44;
	
	if(p->B76==5)
    gcval_press=45;
	
	if(p->B76==0)
    gcval_pcorr=140;

    if(p->B76==1)
    gcval_pcorr=141;

    if(p->B76==2)
    gcval_pcorr=142;

    if(p->B76==3)
    gcval_pcorr=143;
	
	if(p->B76==4)
    gcval_pcorr=144;
	
	if(p->B76==5)
    gcval_pcorr=145;
	
	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
}

void simple::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
	
}
