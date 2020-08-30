/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"momentum_FSI.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

momentum_FSI::momentum_FSI
(
    lexer *p, 
    fdm *a, 
    convection *pconvection, 
    diffusion *pdiffusion, 
    pressure* ppressure, 
    poisson* ppoisson,
    turbulence *pturbulence, 
    solver *psolver, 
    solver *ppoissonsolver, 
    ioflow *pioflow
)
:bcmom(p),un(p),vn(p),wn(p),Fp(p),Fn(p),Fnn(p),Fnnn(p),Gp(p),Gn(p),Gnn(p),Gnnn(p),Hp(p),Hn(p),Hnn(p),Hnnn(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
	ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
}

int momentum_FSI::pcfInd = 0;

momentum_FSI::~momentum_FSI()
{
}

void momentum_FSI::start(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom)
{	
    if (pcfInd == 0) 
    {        
        pflow->discharge(p, a, pgc);
        pflow->inflow(p, a, pgc, a->u, a->v, a->w); 

		//- Velocity predictor step
        predictorStep(p, a, pgc);  
    }
    else if (pcfInd == 1) 
    {
		//- Velocity corrector step
		correctorStep(p,a,pgc);       
    }
	else
    {
		//- Pressure solution
		pflow->pressure_io(p,a,pgc);
		ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,1.0);


		pflow->u_relax(p,a,pgc,a->u);
		pflow->v_relax(p,a,pgc,a->v);
		pflow->w_relax(p,a,pgc,a->w);
		pflow->p_relax(p,a,pgc,a->press);

		cout<<"Call after pressure"<<endl;
		pgc->start1(p,a->u,gcval_u);
		pgc->start2(p,a->v,gcval_v);
		pgc->start3(p,a->w,gcval_w);

		// Save time step
		fieldtimesave(p,a,pgc,pmom); 
    }
}

void momentum_FSI::predictorStep(lexer *p, fdm *a, ghostcell *pgc)
{	
    // Get predictor fluxes
    getF(p, a, pgc, a->u, a->v, a->w);
    getG(p, a, pgc, a->u, a->v, a->w);
    getH(p, a, pgc, a->u, a->v, a->w);
    
	
	// Set fields for first time steps
	if (p->count == 1)
    {
        ULOOP
        {
            Fn(i,j,k) = a->F(i,j,k);
            Fnn(i,j,k) = Fn(i,j,k);
            Fnnn(i,j,k) = Fnn(i,j,k);
		}
		VLOOP
		{
			Gn(i,j,k) = a->G(i,j,k);
			Gnn(i,j,k) = Gn(i,j,k);
			Gnnn(i,j,k) = Gnn(i,j,k);
		}
		WLOOP
		{
			Hn(i,j,k) = a->H(i,j,k);
			Hnn(i,j,k) = Hn(i,j,k);
			Hnnn(i,j,k) = Hnn(i,j,k);
		}
       
        dtn = p->dt;
        dtnn = dtn;
        dtnnn = dtnn;
    }
	
	
    // Coefficients for 4th-order Adam-Bashford
    double beta0 = 
        (3.0*pow(p->dt,3)+4.0*pow(p->dt,2)*(3.0*dtn+2.0*dtnn+dtnnn) + 
        6.0*p->dt*(dtn*(3.0*dtn+4.0*dtnn+2.0*dtnnn)+dtnn*(dtnn+dtnnn)) + 12.0*dtn*(dtn*(dtn+2.0*dtnn+dtnnn)+dtnn*(dtnn+dtnnn)))
        /(12.0*dtn*(dtn+dtnn)*(dtn+dtnn+dtnnn));
    double beta1 = 
        -p->dt*(3.0*pow(p->dt,2)+4.0*p->dt*(2.0*dtn+2.0*dtnn+dtnnn) + 6.0*dtn*(dtn+2.0*dtnn+dtnnn)+6.0*dtnn*(dtnn+dtnnn))
        /(12.0*dtn*dtnn*(dtnn)+dtnnn);
    double beta2 = (3.0*pow(p->dt,2)+4.0*p->dt*(2.0*dtn+dtnn+dtnnn) + 6*dtn*(dtn+dtnn+dtnnn))*p->dt/(12.0*dtnn*dtnnn*(dtn)+dtnn);
    double beta3 = -(3.0*pow(p->dt,2)+4.0*p->dt*(2.0*dtn+dtnn) + 6.0*dtn*(dtn+dtnn))*p->dt/(12.0*dtnnn*(dtnn+dtnnn)*(dtn+dtnn+dtnnn));  
	
	
    // U
    starttime=pgc->timer();
    
	ULOOP
	{
		un(i,j,k) = a->u(i,j,k);
		Fp(i,j,k) = a->F(i,j,k);
		
		a->u(i,j,k) += p->dt*CPOR1*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->F(i,j,k) - (p->dt/p->dt_old)*Fn(i,j,k));

	//	a->u(i,j,k) += p->dt*CPOR1*(a->F(i,j,k)*beta3 + Fn(i,j,k)*beta2 + Fnn(i,j,k)*beta1 + Fnnn(i,j,k)*beta0);
    }
	pgc->start1(p,a->u,gcval_u);
	
	p->utime+=pgc->timer()-starttime;


    // V 
    starttime=pgc->timer();
	VLOOP
	{
		vn(i,j,k) = a->v(i,j,k);
		Gp(i,j,k) = a->G(i,j,k); 
		
		a->v(i,j,k) += p->dt*CPOR2*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->G(i,j,k) - (p->dt/p->dt_old)*Gn(i,j,k));
	
	//	a->v(i,j,k) += p->dt*CPOR1*(a->G(i,j,k)*beta3 + Gn(i,j,k)*beta2 + Gnn(i,j,k)*beta1 + Gnnn(i,j,k)*beta0);
	}
	pgc->start2(p,a->v,gcval_v);
	
	p->vtime+=pgc->timer()-starttime;


    // W 
    starttime=pgc->timer();	
	WLOOP
	{
		wn(i,j,k) = a->w(i,j,k);
		Hp(i,j,k) = a->H(i,j,k);
		
		a->w(i,j,k) += p->dt*CPOR3*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->H(i,j,k) - (p->dt/p->dt_old)*Hn(i,j,k));	
		
	//	a->w(i,j,k) += p->dt*CPOR1*(a->H(i,j,k)*beta3 + Hn(i,j,k)*beta2 + Hnn(i,j,k)*beta1 + Hnnn(i,j,k)*beta0);
    }
	cout<<"Call after predictor"<<endl;
	pgc->start3(p,a->w,gcval_w);
	
	p->wtime+=pgc->timer()-starttime;


// second time --> necessary?
	cout<<"Call after predictor 2"<<endl;
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
}

void momentum_FSI::correctorStep(lexer *p, fdm *a, ghostcell *pgc)
{
    double zeta1 = (3.0*pow(p->dt,2)+4.0*p->dt*(2.0*dtn+dtnn) + 6.0*dtn*(dtn+dtnn))/(12.0*(p->dt+dtn)*(p->dt+dtn+dtnn));
    double zeta2 = (pow(p->dt,2)+2.0*p->dt*(2.0*dtn+dtnn) + 6.0*dtn*(dtn+dtnn))/(12.0*dtn*(dtn+dtnn));
    double zeta3 = (3.0*p->dt-2.0*(2.0*p->dt+dtn+dtnn))*pow(p->dt,2)/(12.0*dtn*dtnn*(p->dt+dtn));
    double zeta4 = (p->dt+2.0*dtn)*pow(p->dt,2)/(12.0*dtnn*(dtn+dtnn)*(p->dt+dtn+dtnn));
     
   
    // U
    starttime=pgc->timer();
   
	for (int qq=0; qq<3; qq++)
    {
        getF(p,a,pgc,a->u,a->v,a->w);
        
        ULOOP
        {
            a->u(i,j,k) = un(i,j,k) + p->dt*(a->F(i,j,k)*zeta4 + Fp(i,j,k)*zeta3 + Fn(i,j,k)*zeta2 + Fnn(i,j,k)*zeta1);
        }    
    }
	pgc->start1(p,a->u,gcval_u);
	
	p->utime+=pgc->timer()-starttime;


    //V 
    starttime=pgc->timer();

	for (int qq=0; qq<3; qq++)
	{
        getG(p,a,pgc,a->u,a->v,a->w);
        
        VLOOP
        {
            a->v(i,j,k) = vn(i,j,k) + p->dt*(a->G(i,j,k)*zeta4 + Gp(i,j,k)*zeta3 + Gn(i,j,k)*zeta2 + Gnn(i,j,k)*zeta1);
        }    
    }
	pgc->start2(p,a->v,gcval_v);
	
	p->vtime+=pgc->timer()-starttime;


    //W 
    starttime=pgc->timer();

	for (int qq=0; qq<3; qq++)
	{
        getH(p,a,pgc,a->u,a->v,a->w);
        
        WLOOP
        {
            a->w(i,j,k) = wn(i,j,k) + p->dt*(a->H(i,j,k)*zeta4 + Hp(i,j,k)*zeta3 + Hn(i,j,k)*zeta2 + Hnn(i,j,k)*zeta1);
        }    
    }
	cout<<"Call after corrector"<<endl;
	pgc->start3(p,a->w,gcval_w);
	
	p->wtime+=pgc->timer()-starttime;


// second time --> necessary?
	cout<<"Call after corrector 2"<<endl;
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);    
}
    
void momentum_FSI::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n]+ a->gi),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FSI::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]+ a->gj),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;	
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FSI::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n]+ a->gk),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}

void momentum_FSI::getF(lexer *p, fdm *a, ghostcell *pgc, field &uvel, field &vvel, field &wvel)
{
    starttime=pgc->timer();

	ULOOP
	{
		a->F(i,j,k) = 0.0;
	}

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,uvel,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,uvel,uvel,vvel,wvel,1.0);
	pconvec->start(p,a,uvel,1,uvel,vvel,wvel);
	pdiff->diff_u(p,a,pgc,psolv,uvel,vvel,wvel,1.0);    							// adds implicit diffusion
    
	p->utime=pgc->timer()-starttime;
} 

void momentum_FSI::getG(lexer *p, fdm *a, ghostcell *pgc, field &uvel, field &vvel, field &wvel)
{
    starttime=pgc->timer();

	VLOOP
	{
		a->G(i,j,k) = 0.0;
	}

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,vvel,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vvel,uvel,vvel,wvel,1.0);
	pconvec->start(p,a,vvel,2,uvel,vvel,wvel);
	pdiff->diff_v(p,a,pgc,psolv,uvel,vvel,wvel,1.0); 

	p->vtime=pgc->timer()-starttime; 
} 

void momentum_FSI::getH(lexer *p, fdm *a, ghostcell *pgc, field &uvel, field &vvel, field &wvel)
{
	starttime=pgc->timer();

	WLOOP
	{
		a->H(i,j,k) = 0.0;
	}

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,wvel,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wvel,uvel,vvel,wvel,1.0);	
	pconvec->start(p,a,wvel,3,uvel,vvel,wvel);	
	pdiff->diff_w(p,a,pgc,psolv,uvel,vvel,wvel,1.0); 
  
 
WLOOP
{
if (fabs(p->pos_x()-0.025)<0.001 && fabs(p->pos_z()-1.05)<0.3) cout<<"Flux "<<p->pos3_x()<<" "<<p->pos3_z()<<" "<<a->H(i,j,k)<<endl;	
}
 
  
	p->wtime=pgc->timer()-starttime;
} 


void momentum_FSI::utimesave(lexer *p, fdm *a, ghostcell* pgc)
{
    ULOOP
    {
		Fnnn(i,j,k) = Fnn(i,j,k);
		Fnn(i,j,k) = Fn(i,j,k);
		Fn(i,j,k) = a->F(i,j,k);
	}
}

void momentum_FSI::vtimesave(lexer *p, fdm *a, ghostcell* pgc)
{
    VLOOP
    {
		Gnnn(i,j,k) = Gnn(i,j,k);
		Gnn(i,j,k) = Gn(i,j,k);
		Gn(i,j,k) = a->G(i,j,k);
    }
}

void momentum_FSI::wtimesave(lexer *p, fdm *a, ghostcell* pgc)
{
    WLOOP
    {
		Hnnn(i,j,k) = Hnn(i,j,k);
		Hnn(i,j,k) = Hn(i,j,k);
		Hn(i,j,k) = a->H(i,j,k);
    }
}

void momentum_FSI::fieldtimesave(lexer *p, fdm *a, ghostcell* pgc, momentum *pmom)
{
    dtnnn = dtnn;
    dtnn = dtn;
    dtn = p->dt;    
}


