/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"pjm_hydrostatic.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
#include"density.h"
 
pjm_hydrostatic::pjm_hydrostatic(density* ppd)
{
    pd = ppd;
    
    gcval_press=40;  
	
	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
}

pjm_hydrostatic::~pjm_hydrostatic()
{
}

void pjm_hydrostatic::start(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, solver* psolv, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);

	pgc->start4(p,a->press,gcval_press);
	
	ucorr(p,a,uvel,alpha);
	vcorr(p,a,vvel,alpha);
	wcorr(p,a,wvel,alpha);
    
    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void pjm_hydrostatic::ucorr(lexer* p, fdm* a, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*pd->roface(p,a,1,0,0)));
}

void pjm_hydrostatic::vcorr(lexer* p, fdm* a, field& vvel,double alpha)
{	
    VLOOP
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))
    /(p->DYP[JP]*(pd->roface(p,a,0,1,0)));
}

void pjm_hydrostatic::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{	
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
}
 
void pjm_hydrostatic::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w, double alpha)
{
    double H,roval,phival,epsi,psi;
    
    
    if(p->j_dir==0)        
    psi = 1.6*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    psi = 1.6*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);

    LOOP
    {
        phival = a->phi(i,j,k);

        if(phival>psi)
            H=1.0;
        else if(phival<-psi)
            H=0.0;
        else
            H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
        
        roval = p->W1*H + p->W3*(1.0-H);
        
        a->press(i,j,k) = a->phi(i,j,k)*roval*fabs(p->W22); 
    } 
}
 
void pjm_hydrostatic::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
}

void pjm_hydrostatic::upgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_hydrostatic::vpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_hydrostatic::wpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_hydrostatic::ini(lexer*p,fdm* a, ghostcell *pgc)
{

}

