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

#include"pjm_IMEX.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"density.h"
#include"ioflow.h"
 
pjm_IMEX::pjm_IMEX(lexer* p, fdm *a) : density(p), pcorr(p), Fp(p), pressn(p)
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
	
	gcval_u=7;
	gcval_v=8;
	gcval_w=9;

    omega_k = 0.0;
}

pjm_IMEX::~pjm_IMEX()
{
}

void pjm_IMEX::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
    
    LOOP
    pcorr(i,j,k)=0.0;
    pgc->start4(p,pcorr,40);
	
    ppois->estart(p,a,pcorr);
	
        starttime=pgc->timer();

    psolv->start(p,a,pgc,pcorr,a->xvec,a->rhsvec,5,gcval_press,p->N44);
	
        endtime=pgc->timer();
    
    pgc->start4(p,pcorr,40);
    presscorr(p,a,pgc,uvel,vvel,wvel,pcorr,alpha);
	pgc->start4(p,a->press,gcval_press);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
   
    pressure_norm(p,a,pgc);
}

void pjm_IMEX::ucorr(fdm* a, lexer* p, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((pcorr(i+1,j,k)-pcorr(i,j,k))
	/(p->DXP[IP]*roface(p,a,1,0,0)));
}

void pjm_IMEX::vcorr(fdm* a, lexer* p, field& vvel,double alpha)
{	
	VLOOP
	vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*((pcorr(i,j+1,k)-pcorr(i,j,k))
	/(p->DYP[JP]*roface(p,a,0,1,0)));
}

void pjm_IMEX::wcorr(fdm* a, lexer* p, field& wvel,double alpha)
{	
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((pcorr(i,j,k+1)-pcorr(i,j,k))
	/(p->DZP[KP]*roface(p,a,0,0,1)));
}

void pjm_IMEX::presscorr(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel, field& pcorr, double alpha)
{
    double rhoU, rhoV, rhoW, sum1, sum2;

    sum1 = 0.0;
    sum2 = 0.0;

    // Calculate Fp as rhs of pressure correction
    LOOP
    {
        Fp(i,j,k) = a->press(i,j,k) + pcorr(i,j,k);
    
        if (alpha == 0.5)
        {
            rhoU = roface(p,a,1,0,0)*uvel(i,j,k);
            i--;
            rhoU -= roface(p,a,1,0,0)*uvel(i,j,k);
            i++;

            rhoV = roface(p,a,0,1,0)*vvel(i,j,k);
            j--;
            rhoV -= roface(p,a,0,1,0)*vvel(i,j,k);
            j++;

            rhoW = roface(p,a,0,0,1)*wvel(i,j,k);
            k--;
            rhoW -= roface(p,a,0,0,1)*wvel(i,j,k);
            k++;

            Fp(i,j,k) -= (a->visc(i,j,k) + a->eddyv(i,j,k))*(rhoU/p->DXN[IP] + rhoV/p->DYN[JP] + rhoW/p->DZN[KP]);
        }
    }
       
    // Calculate relaxation factor omega_k
    LOOP
    {
        sum1 += (pressn(i,j,k) - a->press(i,j,k))*(a->press(i,j,k) - Fp(i,j,k));
        sum2 += ((pressn(i,j,k) - a->press(i,j,k))*(pressn(i,j,k) - 2.0*a->press(i,j,k) + Fp(i,j,k)));
    }

    sum1 = pgc->globalsum(sum1);
    sum2 = pgc->globalsum(sum2);

    omega_k = -sum1/sum2;

    omega_k = 0.0;

    LOOP
    {
        // Store old pressure
        pressn(i,j,k) = a->press(i,j,k);
        
        // Correct pressure
        a->press(i,j,k) = (1.0 - omega_k)*Fp(i,j,k) + omega_k*a->press(i,j,k);
    }
}
 
void pjm_IMEX::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
    count=0;
	double uvel,vvel,wvel;
	
    NLOOP4
	a->rhsvec.V[n]=0.0;
	
    pip=p->Y50;

    LOOP
    {
    a->rhsvec.V[count] =  
                       -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
					   -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])
					   -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
                           
    ++count;
    }
    
    pip=0;
}
 

void pjm_IMEX::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}

void pjm_IMEX::pressure_norm(lexer*p, fdm* a, ghostcell* pgc)
{
    double sum1 = 0.0;
    double sum2 = 0.0;

    LOOP
    {
        sum1 += ((pcorr(i,j,k)*pcorr(i,j,k)));
        sum2 += (a->press(i,j,k)*a->press(i,j,k));
    }

    sum1 = pgc->globalsum(sum1);
    sum2 = pgc->globalsum(sum2);

    cout<<"Pressure splitting error = "<<" "<<sqrt(sum1/sum2)<<endl;
}

void pjm_IMEX::upgrad(lexer*p,fdm* a)
{
}

void pjm_IMEX::vpgrad(lexer*p,fdm* a)
{
}

void pjm_IMEX::wpgrad(lexer*p,fdm* a)
{
}

void pjm_IMEX::fillapu(lexer*p,fdm* a)
{
}

void pjm_IMEX::fillapv(lexer*p,fdm* a)
{
}

void pjm_IMEX::fillapw(lexer*p,fdm* a)
{
}

void pjm_IMEX::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}






