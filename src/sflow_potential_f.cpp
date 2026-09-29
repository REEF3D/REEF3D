/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
    
    // pure Neumann problem: the outflow flux is scaled to the inflow discharge (compatibility)
    Qin2D(p,b,pgc);
    Qout2D(p,b,pgc);
    
    fac_i = 1.0;
    fac_o = fabs(Qo_pf)>1.0e-20?Qi_pf/Qo_pf:0.0;
    
    if(p->mpirank==0)
    cout<<"Qi_pf: "<<Qi_pf<<" Qo_pf: "<<Qo_pf<<" fac_o: "<<fac_o<<endl;

    pgc->gcsl_start4(p,psi,gcval_pot);
    
    laplace(p,b,psi);
    psolv->start(p,pgc,psi,b->M,b->xvec,b->rhsvec,4);
    pgc->gcsl_start4(p,psi,gcval_pot);
    
    p->laplaceiter=p->solveriter;
    
    // cell-centred velocities from the unit discharge q = grad(psi)
    ucalc(p,b,psi);
	vcalc(p,b,psi);

    endtime=pgc->timer();
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0)
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
            b->rhsvec.V[n] += b->M.s[n]*(fac_i*p->Ui*HP)*p->DXP[IM1];
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
            b->rhsvec.V[n] -= b->M.n[n]*(fac_o*p->Uo*HP)*p->DXP[IP1];
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
    // cell-centred U = 0.5*(q_w + q_e)/WL, boundary faces carry the prescribed discharge
    double qw,qe;
    
    SLICELOOP4
    {
    b->U(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        qw = 0.0;
        qe = 0.0;
        
        if(p->flagslice4[Im1J]>0 && p->wet[Im1J]==1)
        qw = (phi(i,j)-phi(i-1,j))/p->DXP[IM1];
        
        else
        if(bc(i-1,j)==1)
        qw = fac_i*p->Ui*HP;
        
        if(p->flagslice4[Ip1J]>0 && p->wet[Ip1J]==1)
        qe = (phi(i+1,j)-phi(i,j))/p->DXP[IP];
        
        else
        if(bc(i+1,j)==2)
        qe = fac_o*p->Uo*HP;
        
        b->U(i,j) = 0.5*(qw+qe)/(b->WL(i,j)>p->A244?b->WL(i,j):1.0e20);
        }
    }
}

void sflow_potential_f::vcalc(lexer *p, fdm2D *b, slice &phi)
{	
    // lateral boundaries are walls
    double qs,qn;
    
    SLICELOOP4
    {
    b->V(i,j) = 0.0;
    
        if(p->wet[IJ]==1 && p->j_dir==1)
        {
        qs = 0.0;
        qn = 0.0;
        
        if(p->flagslice4[IJm1]>0 && p->wet[IJm1]==1)
        qs = (phi(i,j)-phi(i,j-1))/p->DYP[JM1];
        
        if(p->flagslice4[IJp1]>0 && p->wet[IJp1]==1)
        qn = (phi(i,j+1)-phi(i,j))/p->DYP[JP];
        
        b->V(i,j) = 0.5*(qs+qn)/(b->WL(i,j)>p->A244?b->WL(i,j):1.0e20);
        }
    }
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
    // prescribed discharge through the inflow faces (Neumann flux of the Laplace problem)
    Ai=0.0;
    Qi_pf=0.0;

    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        if(p->wet[IJ]==1)
        {
        area = p->DYN[JP]*b->hp(i,j);
        
        Ai+=area;
        Qi_pf+=area*p->Ui;
        }
    }
    
    Ai=pgc->globalsum(Ai);
    Qi_pf=pgc->globalsum(Qi_pf);
}

void sflow_potential_f::Qout2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // discharge through the outflow faces before scaling
    Ao=0.0;
    Qo_pf=0.0;

    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
        
        if(p->wet[IJ]==1)
        {
        area = p->DYN[JP]*b->hp(i,j);
        
        Ao+=area;
        Qo_pf+=area*p->Uo;
        }
    }
    
    Ao=pgc->globalsum(Ao);
    Qo_pf=pgc->globalsum(Qo_pf);
}
