/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"pjm_4th.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
 
pjm_4th::pjm_4th(lexer* p, fdm *a) : density(p)
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
}

pjm_4th::~pjm_4th()
{
}

void pjm_4th::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
	
    ppois->estart(p,a,a->press);
	
        starttime=pgc->timer();

    psolv->start(p,a,pgc,a->press,a->xvec,a->rhsvec,6,gcval_press,p->N44);
	
        endtime=pgc->timer();
    
	pgc->start4(p,a->press,gcval_press);
	
	ucorr(a,p,uvel,alpha);
	vcorr(a,p,vvel,alpha);
	wcorr(a,p,wvel,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
    
}

void pjm_4th::ucorr(fdm* a, lexer* p, field& uvel,double alpha)
{	
	ULOOP
    {
    grad4 = (-a->press(i+2,j,k) + 27.0*a->press(i+1,j,k) - 27.0*a->press(i,j,k) + a->press(i-1,j,k))
          /(-p->XP[IP2] + 27.0*p->XP[IP1] - 27.0*p->XP[IP] + p->XP[IM1]);
       
	uvel(i,j,k) -= (alpha*p->dt*CPOR1*PORVAL1*grad4)/(roface(p,a,1,0,0));
    }
}

void pjm_4th::vcorr(fdm* a, lexer* p, field& vvel,double alpha)
{	
    VLOOP
    {    
    grad4 = (-a->press(i,j+2,k) + 27.0*a->press(i,j+1,k) - 27.0*a->press(i,j,k) + a->press(i,j-1,k))
          /(-p->YP[JP2] + 27.0*p->YP[JP1] - 27.0*p->YP[JP] + p->YP[JM1]);
          
	vvel(i,j,k) -= (alpha*p->dt*CPOR2*PORVAL2*grad4)/(roface(p,a,0,1,0));
    }
}

void pjm_4th::wcorr(fdm* a, lexer* p, field& wvel,double alpha)
{	
	WLOOP
    {    
    grad4 = (-a->press(i,j,k+2) + 27.0*a->press(i,j,k+1) - 27.0*a->press(i,j,k) + a->press(i,j,k-1))
          /(-p->ZP[KP2] + 27.0*p->ZP[KP1] - 27.0*p->ZP[KP] + p->ZP[KM1]);
          
	wvel(i,j,k) -= (alpha*p->dt*CPOR3*PORVAL3*grad4)/(roface(p,a,0,0,1));
    }
}
 
void pjm_4th::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{

    NLOOP4
	a->rhsvec.V[n]=0.0;
	
    pip=p->Y50;

    count=0;
    LOOP
    {
    a->rhsvec.V[count] =  (- (u(i+1,j,k) + 27.0*u(i,j,k) - 27.0*u(i-1,j,k) + u(i-2,j,k))
                           /(-p->XN[IP1] + 27.0*p->XN[IP] - 27.0*p->XN[IM1] + p->XN[IM2])
                           
                          - (v(i,j+1,k) + 27.0*v(i,j,k) - 27.0*v(i,j-1,k) + v(i,j-2,k))
                           /(-p->YN[JP1] + 27.0*p->YN[JP] - 27.0*p->YN[JM1] + p->YN[JM2])
                           
                          - (w(i,j,k+1) + 27.0*w(i,j,k) - 27.0*w(i,j,k-1) + w(i,j,k-2))
                           /(-p->ZN[KP1] + 27.0*p->ZN[KP] - 27.0*p->ZN[KM1] + p->ZN[KM2]))/(alpha*p->dt);
                        
    ++count;
    }
    /*
    count=0;
    LOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
						   -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])
						   -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
    ++count;
    }*/
    
    pip=0;
}
 
void pjm_4th::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}

void pjm_4th::pressure_norm(lexer*p, fdm* a, ghostcell* pgc)
{
    double sum=0.0;

    LOOP
    sum+=a->press(i,j,k);

    sum=pgc->globalsum(sum);

    sum/=double(p->cellnumtot);

    LOOP
    a->press(i,j,k)-=sum;
}

void pjm_4th::upgrad(lexer*p,fdm* a)
{
}

void pjm_4th::vpgrad(lexer*p,fdm* a)
{
}

void pjm_4th::wpgrad(lexer*p,fdm* a)
{
}

void pjm_4th::fillapu(lexer*p,fdm* a)
{
}

void pjm_4th::fillapv(lexer*p,fdm* a)
{
}

void pjm_4th::fillapw(lexer*p,fdm* a)
{
}

void pjm_4th::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}






