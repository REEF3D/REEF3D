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

#include"sflow_potential_f.h"
#include"solver2D.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"lexer.h"
#include<iomanip>

#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HXIPJ (fabs(b->hx(i+1,j))>1.0e-20?b->hx(i+1,j):1.0e20)

#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)
#define HYIJP (fabs(b->hy(i,j+1))>1.0e-20?b->hy(i,j+1):1.0e20)

#define HXP MAX(0.5*(HXIJ + HXIPJ),2.0*p->DXM)
#define HYP MAX(0.5*(HYIJ + HYIJP),2.0*p->DXM)

sflow_potential_f::sflow_potential_f(lexer* p) : bc(p)
{
    gcval_pot=49;
    
    fac_i=fac_o=1.0;
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


    int itermem=p->N46;
    p->N46=2500;
	

    pgc->gcsl_start4(p,psi,gcval_pot);
    
    for(int qn=0; qn<2; ++qn)
    {
    laplace(p,b,psi);
    psolv->start(p,pgc,psi,b->M,b->xvec,b->rhsvec,4);
    pgc->gcsl_start4(p,psi,gcval_pot);
    
    /*
    SLICEBASELOOP
    b->test(i,j) = psi(i,j);
    
    pgc->gcsl_start4(p,b->test,gcval_pot);*/

	
    ucalc(p,b,psi);
	vcalc(p,b,psi);

	pgc->gcsl_start1(p,b->P,10);
	pgc->gcsl_start2(p,b->Q,11);
    
    
    Qin2D(p,b,pgc);
    Qout2D(p,b,pgc);
    }


    endtime=pgc->timer();
    p->laplaceiter=p->solveriter;
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0  && (p->count%p->P12==0))
	cout<<"lapltime: "<<p->laplacetime<<"  lapiter: "<<p->laplaceiter<<endl<<endl;

    p->N46=itermem;
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
        
        b->rhsvec.V[n] = 0.0;
        
        phi(i,j) = 0.0;

    ++n;
    }
    
    
	n=0;
    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        b->M.p[n] =  1.0/(p->DXP[IP]*p->DXN[IP]) + 1.0/(p->DXP[IM1]*p->DXN[IP])
                   + 1.0/(p->DYP[JP]*p->DYN[JP]) + 1.0/(p->DYP[JM1]*p->DYN[JP]);

        b->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP]);
        b->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP]);

        b->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP]);
        b->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP]);
        
        b->rhsvec.V[n] = 0.0;
        }
	++n;
	}
    
    
    n=0;
	SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {

            if((p->flagslice4[Im1J]<0 || p->wet[Im1J]==0) && bc(i-1,j)==0)
            {
            b->M.p[n] += b->M.s[n];
            b->M.s[n] = 0.0;
            }
            
            if((p->flagslice4[Im1J]<0 || p->wet[Im1J]==0) && bc(i-1,j)==1)
            {
            b->rhsvec.V[n] += b->M.s[n]*(0.5*(fac_i+fac_i)*p->Ui*HP)*p->DXP[IM1];
            b->M.p[n] += b->M.s[n];
            b->M.s[n] = 0.0;
            }
            
            if((p->flagslice4[Ip1J]<0 || p->wet[Ip1J]==0) && bc(i+1,j)==0)
            {
            b->M.p[n] += b->M.n[n];
            b->M.n[n] = 0.0;
            }
            
            if((p->flagslice4[Ip1J]<0 || p->wet[Ip1J]==0) && bc(i+1,j)==2)
            {
            b->rhsvec.V[n] -= b->M.n[n]*(0.5*(fac_i+fac_i)*p->Uo*HP)*p->DXP[IP1];
            b->M.p[n] += b->M.n[n];
            b->M.n[n] = 0.0;
            }
            
            if(p->flagslice4[IJm1]<0 || p->wet[IJm1]==0)
            {
            b->M.p[n] += b->M.e[n];
            b->M.e[n] = 0.0;
            }
            
            if(p->flagslice4[IJp1]<0 || p->wet[IJp1]==0)
            {
            b->M.p[n] += b->M.w[n];
            b->M.w[n] = 0.0;
            }
        }
	++n;
	}
}

void sflow_potential_f::ucalc(lexer *p, fdm2D *b, slice &phi)
{	
	SLICELOOP1
    if(p->wet[IJ]==1 && p->wet[Ip1J]==1)
	b->P(i,j) = (phi(i+1,j)-phi(i,j))/(p->DXP[IP]*HXP);
    
    SLICELOOP1
    if(p->wet[IJ]==0 || p->wet[Ip1J]==0)
	b->P(i,j) = 0.0;
}

void sflow_potential_f::vcalc(lexer *p, fdm2D *b, slice &phi)
{	
	SLICELOOP2
    if(p->wet[IJ]==1 && p->wet[IJp1]==1)
	b->Q(i,j) = (phi(i,j+1)-phi(i,j))/(p->DYP[JP]*HYP);
    
    SLICELOOP2
    if(p->wet[IJ]==0 || p->wet[IJp1]==0)
	b->Q(i,j) = 0.0;
}

void sflow_potential_f::ini_bc(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    bc(i,j)=0;
    
    SLICELOOP4
    {
        if(p->flagslice4[Im1J]<0 || p->wet[Im1J]==0)
		bc(i-1,j)=0;
		
		if(p->flagslice4[Ip1J]<0 || p->wet[Ip1J]==0)
		bc(i+1,j)=0;
		
		if(p->flagslice4[IJm1]<0 || p->wet[IJm1]==0)
		bc(i,j-1)=0;
		
		if(p->flagslice4[IJp1]<0 || p->wet[IJp1]==0)
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

void sflow_potential_f::Qin2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    area=0.0;
    Ai=0.0;
    Qi_pf=0.0;
    Ui_pf=0.0;

    // in
    for(n=0;n<p->gcslin_count;n++)
    {
    area=0.0;
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        if(p->wet[IJ]==1)
        {
    
        area = p->DYN[JP]*b->hp(i-1,j);
        
        Ai+=area;
        Qi_pf+=area*b->P(i,j);
        }
    }
    
    Ai=pgc->globalsum(Ai);
    Qi_pf=pgc->globalsum(Qi_pf);
    
    fac_i = p->Qi/Qi_pf;
 
    if(p->mpirank==0)
    cout<<"Qi_pf: "<<Qi_pf<<" fac_i: "<<fac_i<<endl;
}


void sflow_potential_f::Qout2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    area=0.0;
    Ao=0.0;
    Qo_pf=0.0;
    Uo_pf=0.0;

    // out
    for(n=0;n<p->gcslout_count;n++)
    {
    area=0.0;
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
        
        if(p->wet[IJ]==1)
        {
    
        area = p->DYN[JP]*b->hp(i,j);
        
        Ao+=area;
        Qo_pf+=area*b->P(i,j);
        }
    }
    
    Ao=pgc->globalsum(Ao);
    Qo_pf=pgc->globalsum(Qo_pf);

    fac_o = p->Qi/Qo_pf;
    
    if(p->mpirank==0)
    cout<<"Qo_pf: "<<Qo_pf<<" fac_o: "<<fac_o<<endl;
}