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

#include"sflow_potential_f.h"
#include"solver2D.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"lexer.h"
#include<iomanip>

#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

sflow_potential_f::sflow_potential_f(lexer* p) : bc(p)
{
    gcval_pot=49;
}

sflow_potential_f::~sflow_potential_f()
{
}

void sflow_potential_f::start(lexer *p, fdm2D *b, solver2D *psolv, ghostcell *pgc)
{
    if(p->mpirank==0 )
	cout<<"starting potential flow solver..."<<endl<<endl;
    
    slice4 psi(p);
    
    ini_bc(p,b,pgc);

    starttime=pgc->timer();
	
	pgc->gcsl_start4(p,psi,gcval_pot);
	
	SLICELOOP4
	psi(i,j) = (p->pos_x()-p->global_xmin)*p->Ui;
	
    ucalc(p,b,psi);
	vcalc(p,b,psi);
    
    pgc->gcsl_start1(p,b->P,10);
	pgc->gcsl_start2(p,b->Q,11);


    int itermem=p->N46;
    p->N46=2500;
	

    pgc->gcsl_start4(p,psi,gcval_pot);


    laplace(p,b,psi);
    psolv->start(p,pgc,psi,b->M,b->xvec,b->rhsvec,4,gcval_pot,p->N44);
    pgc->gcsl_start4(p,psi,gcval_pot);

	
    ucalc(p,b,psi);
	vcalc(p,b,psi);


	pgc->gcsl_start1(p,b->P,10);
	pgc->gcsl_start2(p,b->Q,11);


    endtime=pgc->timer();
    p->laplaceiter=p->solveriter;
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0  && (p->count%p->P12==0))
	cout<<"lapltime: "<<p->laplacetime<<"  lapiter: "<<p->laplaceiter<<endl<<endl;

    p->N46=itermem;
}

void sflow_potential_f::ucalc(lexer *p, fdm2D *b, slice &phi)
{	
	SLICELOOP1
	b->P(i,j) = (phi(i+1,j)-phi(i,j))/(p->DXP[IP]*HP);
	

}

void sflow_potential_f::vcalc(lexer *p, fdm2D *b, slice &phi)
{	
	SLICELOOP2
	b->Q(i,j) = (phi(i,j+1)-phi(i,j))/(p->DYP[JP]*HP);


}

void sflow_potential_f::rhs(lexer *p, fdm2D *b)
{
    count=0;
    SLICELOOP4
    {
    b->rhsvec.V[count] = 0.0;
    count++;
    }
}

void sflow_potential_f::laplace(lexer *p, fdm2D *b, slice &phi)
{
    n=0;
    SLICEBASELOOP
    {
    b->M.p[n]  =  1.0;


        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;

        b->M.w[n] = 0.0;
        b->M.e[n] = 0.0;
        
        b->rhsvec.V[n] =  0.0;
    ++n;
    }
    
    
	n=0;
    SLICELOOP4
	{
	b->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP]) + 1.0/(p->DXP[IM1]*p->DXN[IP])
                + 1.0/(p->DYP[JP]*p->DYN[JP]) + 1.0/(p->DYP[JM1]*p->DYN[JP]);

   	b->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP]);
	b->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP]);

	b->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP]);
	b->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP]);
	
	b->rhsvec.V[n] = 0.0;
	
	++n;
	}
    
    
    n=0;
	SLICELOOP4
	{
		if(p->flagslice4[Im1J]<0 && bc(i-1,j)==0)
		{
		b->M.p[n] += b->M.s[n];
		b->M.s[n] = 0.0;
		}
        
        if(p->flagslice4[Im1J]<0 && bc(i-1,j)==1)
		{
        b->rhsvec.V[n] += b->M.s[n]*(p->Ui/HP)*p->DXP[IM1];
		b->M.p[n] += b->M.s[n];
		b->M.s[n] = 0.0;
		}
		
		if(p->flagslice4[Ip1J]<0 && bc(i+1,j)==0)
		{
		b->M.p[n] += b->M.n[n];
		b->M.n[n] = 0.0;
		}
        
        if(p->flagslice4[Ip1J]<0 && bc(i+1,j)==2)
		{
        b->rhsvec.V[n] -= b->M.n[n]*(p->Uo/HP)*p->DXP[IP1];

		b->M.n[n] = 0.0;
		}
		
		if(p->flagslice4[IJm1]<0)
		{
		b->M.p[n] += b->M.e[n];
		b->M.e[n] = 0.0;
		}
		
		if(p->flagslice4[IJp1]<0)
		{
		b->M.p[n] += b->M.w[n];
		b->M.w[n] = 0.0;
		}

	++n;
	}
}

void sflow_potential_f::ini_bc(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    bc(i,j)=0;
    
    SLICELOOP4
    {
        if(p->flagslice4[Im1J]<0)
		bc(i-1,j)=0;
		
		if(p->flagslice4[Ip1J]<0)
		bc(i+1,j)=0;
		
		if(p->flagslice4[IJm1]<0)
		bc(i,j-1)=0;
		
		if(p->flagslice4[IJp1]<0)
		bc(i,j+1)=0;
    }
    

    GCSL4LOOP
    {
    
        if(p->gcbsl4[n][4]==1)
        {
            i = p->gcbsl4[n][0];
            j = p->gcbsl4[n][1];
       
            if(p->gcbsl4[n][3]==1)
            bc(i-1,j)=1;
            
            if(p->gcbsl4[n][3]==3)
            bc(i,j-1)=1;
            
            if(p->gcbsl4[n][3]==2)
            bc(i,j+1)=1;
            
            if(p->gcbsl4[n][3]==4)
            bc(i+1,j)=1;
            
            cout<<"BC1: "<<endl;
 
        }
        
        if(p->gcbsl4[n][4]==2)
        {
            i=p->gcbsl4[n][0]; 
            j=p->gcbsl4[n][1];
      
       
            if(p->gcbsl4[n][3]==1)
            bc(i-1,j)=2;
            
            if(p->gcbsl4[n][3]==3)
            bc(i,j-1)=2;
            
            if(p->gcbsl4[n][3]==2)
            bc(i,j+1)=2;
            
            if(p->gcbsl4[n][3]==4)
            bc(i+1,j)=2;
 
        }
    }
    
}
