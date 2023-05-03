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

#include"potential_water.h"
#include"solver.h"
#include"ghostcell.h"
#include"fdm.h"
#include"lexer.h"
#include<iomanip>

potential_water::potential_water(lexer* p) : bc(p), eps(-1.6*p->DXM)
{
    gcval_pot=49;
}

potential_water::~potential_water()
{
}

void potential_water::start(lexer*p,fdm* a,solver* psolv, ghostcell* pgc)
{
    if(p->mpirank==0 )
	cout<<"starting potential_water flow solver..."<<endl<<endl;
    
    field4 psi(p);
    
    ini_bc(p,a,pgc);

    starttime=pgc->timer();
	
    int itermem=p->N46;
    p->N46=2500;

    laplace(p,a,psi);
	psolv->start(p,a,pgc,psi,a->rhsvec,4);
    pgc->start4(p,psi,gcval_pot);

	
    ucalc(p,a,psi);
	vcalc(p,a,psi);
	wcalc(p,a,psi);

	pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);

    endtime=pgc->timer();
    p->laplaceiter=p->solveriter;
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0  && (p->count%p->P12==0))
	cout<<"lapltime: "<<p->laplacetime<<"  lapiter: "<<p->laplaceiter<<endl<<endl;

    p->N46=itermem;
    
    
    LOOP
    psi(i,j,k) = 0.0;
}

void potential_water::ucalc(lexer *p, fdm *a, field &phi)
{	
	ULOOP
    if(a->phi(i,j,k)>eps && a->phi(i+1,j,k)>eps)
	a->u(i,j,k) = (phi(i+1,j,k)-phi(i,j,k))/p->DXP[IP];
	
	if(p->I21==1)
	ULOOP
	if(0.5*(a->phi(i,j,k)+a->phi(i+1,j,k))<-p->F45*p->DXP[IP])
	a->u(i,j,k)=0.0;
    
    if(p->S10==2)
	ULOOP
	if(0.5*(a->topo(i,j,k)+a->topo(i+1,j,k))<-p->F45*p->DXP[IP])
	a->u(i,j,k)=0.0;
}

void potential_water::vcalc(lexer *p, fdm *a, field &phi)
{	
	VLOOP
    if(a->phi(i,j,k)>eps && a->phi(i,j+1,k)>eps)
	a->v(i,j,k) = (phi(i,j+1,k)-phi(i,j,k))/p->DYP[JP];

	if(p->I21==1)
	VLOOP
	if(a->phi(i,j,k)<-p->F45*p->DYP[JP])
	a->v(i,j,k)=0.0;
    
    if(p->S10==2)
	VLOOP
	if(0.5*(a->topo(i,j,k)+a->topo(i,j+1,k))<-p->F45*p->DYP[JP])
	a->v(i,j,k)=0.0;
}

void potential_water::wcalc(lexer *p, fdm *a, field &phi)
{
	WLOOP
    if(a->phi(i,j,k)>eps && a->phi(i,j,k+1)>eps)
    a->w(i,j,k) = (phi(i,j,k+1)-phi(i,j,k))/p->DZP[KP];
	
    if(p->I21==1)
	WLOOP
	if(a->phi(i,j,k)<-p->F45*p->DZP[KP])
	a->w(i,j,k)=0.0;
    
    if(p->S10==2)
	WLOOP
	if(0.5*(a->topo(i,j,k)+a->topo(i,j,k+1))<-p->F45*p->DZP[KP])
	a->w(i,j,k)=0.0;
}

void potential_water::rhs(lexer *p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count] = 0.0;
    count++;
    }
}

void potential_water::laplace(lexer *p, fdm *a, field &phi)
{
    n=0;
    BASELOOP
    {
    a->M.p[n]  =  1.0;


        a->M.n[n] = 0.0;
        a->M.s[n] = 0.0;

        a->M.w[n] = 0.0;
        a->M.e[n] = 0.0;

        a->M.t[n] = 0.0;
        a->M.b[n] = 0.0;
        
        a->rhsvec.V[n] =  0.0;
    ++n;
    }
    
    
	n=0;
    LOOP
    {
        if(a->phi(i,j,k)>=eps)
        {
        a->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP]) + 1.0/(p->DXP[IM1]*p->DXN[IP])
                    + 1.0/(p->DYP[JP]*p->DYN[JP]) + 1.0/(p->DYP[JM1]*p->DYN[JP])
                    + 1.0/(p->DZP[KP]*p->DZN[KP]) + 1.0/(p->DZP[KM1]*p->DZN[KP]);

        a->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP]);
        a->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP]);

        a->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP]);
        a->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP]);

        a->M.t[n] = -1.0/(p->DZP[KP]*p->DZN[KP]);
        a->M.b[n] = -1.0/(p->DZP[KM1]*p->DZN[KP]);
        
        a->rhsvec.V[n] = 0.0;
        }
	
	++n;
	}
    
    
    n=0;
	LOOP
    {
    if(a->phi(i,j,k)>=eps)
	{
		if((p->flag4[Im1JK]<0 || a->phi(i-1,j,k)<eps) && bc(i-1,j,k)==0)
		{
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
        
        if((p->flag4[Im1JK]<0 || a->phi(i-1,j,k)<eps) && bc(i-1,j,k)==1)
		{
        a->rhsvec.V[n] += a->M.s[n]*p->Ui*p->DXP[IM1];// - a->M.s[n]*phi(i,j,k);
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
		
		if((p->flag4[Ip1JK]<0 || a->phi(i+1,j,k)<eps) && bc(i+1,j,k)==0)
		{
		a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
        
        if((p->flag4[Ip1JK]<0 || a->phi(i+1,j,k)<eps) && bc(i+1,j,k)==2)
		{
        a->rhsvec.V[n] -= a->M.n[n]*p->Uo*p->DXP[IP1];
		//a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0 || a->phi(i,j-1,k)<eps)
		{
		a->M.p[n] += a->M.e[n];
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0 || a->phi(i,j+1,k)<eps)
		{
		a->M.p[n] += a->M.w[n];
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0 || a->phi(i,j,k-1)<eps)
		{
		a->M.p[n] += a->M.b[n];
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0 || a->phi(i,j,k+1)<eps)
		{
		a->M.p[n] += a->M.t[n];
		a->M.t[n] = 0.0;
		}
        }
	++n;
	}
}

void potential_water::ini_bc(lexer *p, fdm *a, ghostcell *pgc)
{
    LOOP
    bc(i,j,k)=0;
    
    LOOP
    {
        if(p->flag4[Im1JK]<0)
		bc(i-1,j,k)=0;
		
		if(p->flag4[Ip1JK]<0)
		bc(i+1,j,k)=0;
		
		if(p->flag4[IJm1K]<0)
		bc(i,j-1,k)=0;
		
		if(p->flag4[IJp1K]<0)
		bc(i,j+1,k)=0;
		
		if(p->flag4[IJKm1]<0)
		bc(i,j,k-1)=0;
		
		if(p->flag4[IJKp1]<0)
		bc(i,j,k+1)=0;
    }
    

    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
            i=p->gcb4[n][0]; 
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];  
       
            if(p->gcb4[n][3]==1)
            bc(i-1,j,k)=1;
            
            if(p->gcb4[n][3]==3)
            bc(i,j-1,k)=1;
            
            if(p->gcb4[n][3]==2)
            bc(i,j+1,k)=1;
            
            if(p->gcb4[n][3]==4)
            bc(i+1,j,k)=1;
 
        }
        
        if(p->gcb4[n][4]==2)
        {
            i=p->gcb4[n][0]; 
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];  
       
            if(p->gcb4[n][3]==1)
            bc(i-1,j,k)=2;
            
            if(p->gcb4[n][3]==3)
            bc(i,j-1,k)=2;
            
            if(p->gcb4[n][3]==2)
            bc(i,j+1,k)=2;
            
            if(p->gcb4[n][3]==4)
            bc(i+1,j,k)=2;
 
        }
    }
    
}
