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

#include"pjm_fsi.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
 
pjm_fsi::pjm_fsi(lexer* p, fdm *a) : density(p)
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

pjm_fsi::~pjm_fsi()
{
}

void pjm_fsi::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0)
    cout<<".";

	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
	
    ppois->estart(p,a,a->press);
	
        starttime=pgc->timer();

    psolv->start(p,a,pgc,a->press,a->xvec,a->rhsvec,4,gcval_press,p->N44);
	
        endtime=pgc->timer();


	pgc->start4(p,a->press,gcval_press);
	
	ucorr(a,p,uvel,alpha);
	vcorr(a,p,vvel,alpha);
	wcorr(a,p,wvel,alpha);

LOOP
{
if (fabs(p->pos_x()-0.025)<0.001 && fabs(p->pos_z()-1.05)<0.3) cout<<"dw/dz corr "<<p->pos_z()<<" "<<a->w(i,j,k)<<endl;	
if (fabs(p->pos_x()-0.025)<0.001 && fabs(p->pos_z()-1.05)<0.3) cout<<"with p "<<p->pos_z()<<" "<<a->press(i,j,k+1)<<" "<<a->press(i,j,k)<<endl;	
}

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && innercounter==p->N50-1)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
    
}

void pjm_fsi::ucorr(fdm* a, lexer* p, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*roface(p,a,1,0,0)));
}

void pjm_fsi::vcorr(fdm* a, lexer* p, field& vvel,double alpha)
{	
	VLOOP
	vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*((a->press(i,j+1,k)-a->press(i,j,k))
	/(p->DYP[JP]*roface(p,a,0,1,0)));
}

void pjm_fsi::wcorr(fdm* a, lexer* p, field& wvel,double alpha)
{	
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*roface(p,a,0,0,1)));
}
 
void pjm_fsi::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	
	double uvel,vvel,wvel;
	
    NLOOP4
	a->rhsvec.V[n]=0.0;
	
    pip=p->Y50;
    
    count=0;
    LOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
						   -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])
						   -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
    ++count;
	
	a->test(i,j,k) = (w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
	
if (fabs(p->pos_x()-0.025)<0.001 && fabs(p->pos_z()-1.075)<0.3) cout<<"dw/dz "<<p->pos_z()<<" "<<w(i,j,k)<<endl;	
	
    }

    pip=0;
}
 

void pjm_fsi::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
   cout<<"Call before pressure solution"<<endl;
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}

void pjm_fsi::pressure_norm(lexer*p, fdm* a, ghostcell* pgc)
{
    double sum=0.0;

    LOOP
    sum+=a->press(i,j,k);

    sum=pgc->globalsum(sum);

    sum/=double(p->cellnumtot);

    LOOP
    a->press(i,j,k)-=sum;
}

void pjm_fsi::upgrad(lexer*p,fdm* a)
{
    ULOOP
    a->F(i,j,k)-=PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*roface(p,a,1,0,0));
}

void pjm_fsi::vpgrad(lexer*p,fdm* a)
{
    VLOOP
    a->G(i,j,k)-=PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))/(p->DYP[JP]*roface(p,a,0,1,0));
}

void pjm_fsi::wpgrad(lexer*p,fdm* a)
{
    WLOOP
    a->H(i,j,k)-=PORVAL3*(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*roface(p,a,0,0,1));
}

void pjm_fsi::fillapu(lexer*p,fdm* a)
{
}

void pjm_fsi::fillapv(lexer*p,fdm* a)
{
}

void pjm_fsi::fillapw(lexer*p,fdm* a)
{
}

void pjm_fsi::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}






