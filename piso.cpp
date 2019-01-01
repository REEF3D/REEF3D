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

#include"piso.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"

piso::piso(lexer* p,fdm *a) : density(p),apu(p),uc(p),apv(p),vc(p),apw(p),wc(p),pcorr(p)
{
    ini_bcval(p,a);
}

piso::~piso()
{
}

void piso::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    pgc->start1(p,apu,14);
    pgc->start2(p,apv,15);
    pgc->start3(p,apw,16);

    LOOP
    pcorr(i,j,k)=0.0;
    pgc->start4(p,pcorr,gcval_pcorr);

    //1st step
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);
	rhs(p,a,pgc,a->u,a->v,a->w,alpha);
    ppois->istart(p,a,apu,apv,apw,pcorr);
        starttime=pgc->timer();
	psolv->start(p,a,pgc,pcorr,a->xvec,a->rhsvec,5,gcval_press,p->N44);
        endtime=pgc->timer();
    p->poissoniter=p->solveriter;
    presscorr(p,a,pcorr);
	pgc->start4(p,a->press,gcval_press);
	pgc->start4(p,pcorr,gcval_pcorr);

	ucorr(a,p,pcorr);
	vcorr(a,p,pcorr);
	wcorr(a,p,pcorr);

	convergence(p,a,pgc);

	pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);

    p->poissoniter+=p->solveriter;
	p->poissontime+=endtime-starttime;

	if(p->mpirank==0 && innercounter==p->N50-1)
	cout<<"pc: "<<p->pcnorm<<" uc: "<<p->ucnorm<<" vc: "<<p->vcnorm<<" wc: "<<p->wcnorm<<endl;


	//2nd step
	for(pisoiter=0;pisoiter<p->D31;++pisoiter)
	{
    LOOP
    pcorr(i,j,k)=0.0;
    pgc->start4(p,pcorr,gcval_pcorr);

    pmom->fillaij1(p,a,pgc,psolv);
    psolv->fillxvec1(p,a,uc);
    //psolv->gcpara_update(p,a->xvec,pgc);
    filluct(p,a);

    pmom->fillaij2(p,a,pgc,psolv);
    psolv->fillxvec2(p,a,vc);
    //psolv->gcpara_update(p,a->xvec,pgc);
    fillvct(p,a);

    pmom->fillaij3(p,a,pgc,psolv);
    psolv->fillxvec3(p,a,wc);
    //psolv->gcpara_update(p,a->xvec,pgc);
    fillwct(p,a);

    pgc->start1(p,uc,17);
	pgc->start2(p,vc,18);
	pgc->start3(p,wc,19);

	rhs(p,a,pgc,uc,vc,wc,alpha);
	ppois->istart(p,a,apu,apv,apw,pcorr);
        starttime2=pgc->timer();
    psolv->start(p,a,pgc,pcorr,a->xvec,a->rhsvec,5,gcval_press,p->N44);
        endtime2=pgc->timer();
    presscorr(p,a,pcorr);
    pgc->start4(p,a->press,gcval_press);
	pgc->start4(p,pcorr,gcval_pcorr);

	ucorr2(a,p,pcorr);
	vcorr2(a,p,pcorr);
	wcorr2(a,p,pcorr);

	convergence(p,a,pgc);

    p->poissoniter+=p->solveriter;
	p->poissontime+=endtime2-starttime2;

	if(p->mpirank==0&& innercounter==p->N50-1)
	{
	cout<<"pc";
	for(cc=0;cc<=pisoiter;++cc)
	cout<<"c";
	cout<<": "<<p->pcnorm;

	cout<<" uc";
	for(cc=0;cc<=pisoiter;++cc)
	cout<<"c";
	cout<<": "<<p->ucnorm;

	cout<<" vc";
	for(cc=0;cc<=pisoiter;++cc)
	cout<<"c";
	cout<<": "<<p->vcnorm;

	cout<<" wc";
	for(cc=0;cc<=pisoiter;++cc)
	cout<<"c";
	cout<<": "<<p->wcnorm<<endl;
	}
	}

	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	{
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
	p->poissoniter=0;
	p->poissontime=0.0;
	}

}

void piso::ucorr(fdm* a, lexer* p, field& pcorr)
{
    p->ucnorm=0.0;

	ULOOP
	{
	a->u(i,j,k) +=uc(i,j,k)=ucc=- (pcorr(i+1,j,k)-pcorr(i,j,k))
	/(p->DXP[IP]*apu(i,j,k)*roface(p,a,1,0,0));
	p->ucnorm+=ucc*ucc;
	}

}

void piso::vcorr(fdm* a, lexer* p, field& pcorr)
{
    p->vcnorm=0.0;

	VLOOP
	{
	a->v(i,j,k) +=vc(i,j,k)=vcc= -(pcorr(i,j+1,k)-pcorr(i,j,k))
	/(p->DYP[JP]*apv(i,j,k)*roface(p,a,0,1,0));
	p->vcnorm+=vcc*vcc;
	}
}

void piso::wcorr(fdm* a, lexer* p, field& pcorr)
{
    p->wcnorm=0.0;

	WLOOP
	{
	a->w(i,j,k) +=wc(i,j,k)=wcc=- (pcorr(i,j,k+1)-pcorr(i,j,k))
	/(p->DZP[KP]*apw(i,j,k)*roface(p,a,0,0,1));
	p->wcnorm+=wcc*wcc;
	}
}

void piso::ucorr2(fdm* a, lexer* p, field& pcorr)
{
    p->ucnorm=0.0;

	ULOOP
	{
	a->u(i,j,k) +=uc(i,j,k)=ucc= uc(i,j,k)-(pcorr(i+1,j,k)-pcorr(i,j,k))	/(p->DXP[IP]*apu(i,j,k)*0.5*(a->ro(i+1,j,k)+a->ro(i,j,k)));
	p->ucnorm+=ucc*ucc;
	}
}

void piso::vcorr2(fdm* a, lexer* p, field& pcorr)
{
    p->vcnorm=0.0;

	VLOOP
	{
	a->v(i,j,k) +=vc(i,j,k)=vcc= vc(i,j,k)-(pcorr(i,j+1,k)-pcorr(i,j,k))	/(p->DYP[JP]*apv(i,j,k)*0.5*(a->ro(i,j+1,k)+a->ro(i,j,k)));
	p->vcnorm+=vcc*vcc;
	}
}

void piso::wcorr2(fdm* a, lexer* p, field& pcorr)
{
    p->wcnorm=0.0;

	WLOOP
	{
	a->w(i,j,k) +=wc(i,j,k)=wcc= wc(i,j,k)-(pcorr(i,j,k+1)-pcorr(i,j,k)) /(p->DZP[KP]*apw(i,j,k)*0.5*(a->ro(i,j,k+1)+a->ro(i,j,k)));
	p->wcnorm+=wcc*wcc;
	}
}

void piso::presscorr(lexer* p, fdm* a, field& pcorr)
{
    p->pcnorm=0.0;
    LOOP
    {
    a->press(i,j,k)+=pcc=p->N53*pcorr(i,j,k);
    p->pcnorm+=pcc*pcc;
    }
}

void piso::rhs(lexer *p, fdm* a, ghostcell *pgc, field& u, field& v, field& w, double alpha)
{
    pip=p->Y50;

    LOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
						   -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])
						   -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
                           
    ++count;
    }
    
    pip=0;
}

void piso::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}


void piso::upgrad(lexer*p,fdm* a)
{	
    ULOOP
	a->F(i,j,k)-=PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*0.5*(a->ro(i+1,j,k)+a->ro(i,j,k)));
}

void piso::vpgrad(lexer*p,fdm* a)
{	
    VLOOP
	a->G(i,j,k)-=PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))/(p->DYP[JP]*0.5*(a->ro(i,j+1,k)+a->ro(i,j,k)));
}

void piso::wpgrad(lexer*p,fdm* a)
{	
    WLOOP
	a->H(i,j,k)-=PORVAL3*(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*0.5*(a->ro(i,j,k+1)+a->ro(i,j,k)));
}

void piso::fillapu(lexer*p,fdm* a)
{
    count=0;

    ULOOP
    {
        apu(i,j,k)=a->M.p[count];

	++count;
    }
}


void piso::fillapv(lexer*p,fdm* a)
{
    count=0;
	
    VLOOP
    {
        apv(i,j,k)=a->M.p[count];


	++count;
    }
}

void piso::fillapw(lexer*p,fdm* a)
{
    count=0;
	
    WLOOP
    {
        apw(i,j,k)=a->M.p[count];

	++count;
    }
}

void piso::filluct(lexer*p,fdm* a)
{
    n=0;
	
    ULOOP
    {
		
		uc(i,j,k)= a->M.s[n]*a->xvec.V[Im1_J_K_1]
				 + a->M.w[n]*a->xvec.V[I_Jp1_K_1]
				 + a->M.e[n]*a->xvec.V[I_Jm1_K_1]
				 + a->M.t[n]*a->xvec.V[I_J_Kp1_1]
				 + a->M.b[n]*a->xvec.V[I_J_Km1_1];
	++n;
    }

}

void piso::fillvct(lexer*p,fdm* a)
{
    n=0;
	
    VLOOP
    {
        vc(i,j,k)= a->M.s[n]*a->xvec.V[Im1_J_K_2]
				 + a->M.w[n]*a->xvec.V[I_Jp1_K_2]
				 + a->M.e[n]*a->xvec.V[I_Jm1_K_2]
				 + a->M.t[n]*a->xvec.V[I_J_Kp1_2]
				 + a->M.b[n]*a->xvec.V[I_J_Km1_2];

	++n;
    }
}

void piso::fillwct(lexer*p,fdm* a)
{
    n=0;
	
    WLOOP
    {
        wc(i,j,k)= a->M.s[n]*a->xvec.V[Im1_J_K_3]
				 + a->M.w[n]*a->xvec.V[I_Jp1_K_3]
				 + a->M.e[n]*a->xvec.V[I_Jm1_K_3]
				 + a->M.t[n]*a->xvec.V[I_J_Kp1_3]
				 + a->M.b[n]*a->xvec.V[I_J_Km1_3];

	++n;
    }
}

void piso::convergence(lexer *p, fdm *a, ghostcell *pgc)
{
    p->pcnorm = sqrt(pgc->globalsum(p->pcnorm)/double(p->cellnumtot))/(fabs(p->pressmax-p->pressmin)>0.0?fabs(p->pressmax-p->pressmin):1.0);

    p->ucnorm = sqrt(pgc->globalsum(p->ucnorm)/double(p->cellnumtot));
	p->vcnorm = sqrt(pgc->globalsum(p->vcnorm)/double(p->cellnumtot));
	p->wcnorm = sqrt(pgc->globalsum(p->wcnorm)/double(p->cellnumtot));

}

void piso::ini_bcval(lexer *p, fdm *a)
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

void piso::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

