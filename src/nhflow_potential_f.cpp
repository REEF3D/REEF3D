/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_potential_f.h"
#include"solver.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include<iomanip>

nhflow_potential_f::nhflow_potential_f(lexer* p) 
{
    gcval_pot=49;
}

nhflow_potential_f::~nhflow_potential_f()
{
}

void nhflow_potential_f::start(lexer*p, fdm_nhf *d, solver* psolv, ghostcell* pgc)
{
    int itermem;
    p->Darray(PSI,p->imax*p->jmax*(p->kmax+2));
    
    if(p->mpirank==0 )
	cout<<"potential flow solver..."<<endl<<endl;
    
    if(fabs(p->Ui)<1.0e-9)
    {
    if(p->mpirank==0)
    cout<<"no inflow...potential flow stop"<<endl;
    goto finalize;
    }
    
    if(fabs(p->Uo)<1.0e-9)
    {
    if(p->mpirank==0)
    cout<<"no outflow...potential flow stop"<<endl;
    goto finalize;
    }


    starttime=pgc->timer();
	
	pgc->start49V(p,PSI,gcval_pot);
	
	LOOP
	PSI[IJK] = 0.0;
    
    pgc->start49V(p,PSI,gcval_pot);
	
    itermem=p->N46;
    p->N46=5000;
	
    for(int qn=0; qn<1;++qn)
    {
    laplace(p,d,pgc);
    psolv->startV(p,pgc,PSI,d->rhsvec,d->M,44);
    }
    
    pgc->start49V(p,PSI,gcval_pot);
    
    ucalc(p,d);
	vcalc(p,d);
	wcalc(p,d);

	pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);

    endtime=pgc->timer();
    p->laplaceiter=p->solveriter;
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0  && (p->count%p->P12==0))
	cout<<"lapltime: "<<p->laplacetime<<"  lapiter: "<<p->laplaceiter<<endl<<endl;

    p->N46=itermem;
    

    finalize:
    
    p->del_Darray(PSI,p->imax*p->jmax*(p->kmax+2));
}

void nhflow_potential_f::ucalc(lexer *p, fdm_nhf *d)
{	
	LOOP
    {
	d->U[IJK] =  (PSI[Ip1JK]-PSI[Im1JK])/(p->DXP[IP]+p->DXP[IM1]) 
    
              + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KM1]));
              
    d->UH[IJK] = d->WL(i,j)*d->U[IJK];
    }
    
    
    /*
    if(p->X10==1)
	ULOOP
    if(a->fb(i+1,j,k)<0.0 || a->fb(i,j,k)<0.0)
	a->u(i,j,k)=0.0;
    
	ULOOP
    if(p->flagsf4[Ip1JK]<0 || p->flagsf4[IJK]<0.0)
	a->u(i,j,k)=0.0;*/
}

void nhflow_potential_f::vcalc(lexer *p, fdm_nhf *d)
{	
    LOOP
    {
	d->V[IJK] =  (PSI[IJp1K]-PSI[IJm1K])/(p->DYP[JP]+p->DYP[JP1]) 
    
              + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KM1]));
              
    d->VH[IJK] = d->WL(i,j)*d->V[IJK];
    }

    /*
    if(p->X10==1)
	VLOOP
    if(a->fb(i,j+1,k)<0.0 || a->fb(i,j,k)<0.0)
	a->v(i,j,k)=0.0;

	VLOOP
    if(p->flagsf4[IJp1K]<0 || p->flagsf4[IJK]<0.0)
	a->v(i,j,k)=0.0;*/
}

void nhflow_potential_f::wcalc(lexer *p, fdm_nhf *d)
{
    LOOP
    {
	d->W[IJK] =  p->sigz[IJ]*(PSI[IJKp1]-PSI[IJKp1])/(p->DZP[KP]+p->DZP[KP1]); 
              
    d->WH[IJK] = d->WL(i,j)*d->W[IJK];
    }
    
    /*
    if(p->X10==1)
	WLOOP
    if(a->fb(i,j,k+1)<0.0 || a->fb(i,j,k)<0.0)
	a->w(i,j,k)=0.0;
    
	WLOOP
    if(p->flagsf4[IJKp1]<0 || p->flagsf4[IJK]<0.0)
	a->w(i,j,k)=0.0;
    */
}

void nhflow_potential_f::rhs(lexer *p, fdm_nhf *d)
{/*
    count=0;
    LOOP
    {
    a->rhsvec.V[count] = 0.0;
    count++;
    }*/
}

void nhflow_potential_f::laplace(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    n=0;
    BASELOOP
    {
    d->M.p[n]  =  1.0;

        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

        d->M.w[n] = 0.0;
        d->M.e[n] = 0.0;
        
        d->M.t[n] = 0.0;
        d->M.b[n] = 0.0;
        
        d->rhsvec.V[n] =  0.0;
        
    ++n;
    }
    
	
    
    
    n=0;
    LOOP
	{
        WETDRYDEEP
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            
            d->M.p[n]  =  1.0/(p->W1*p->DXP[IP]*p->DXN[IP])
                        + 1.0/(p->W1*p->DXP[IM1]*p->DXN[IP])
                        
                        + 1.0/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + 1.0/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + sigxyz2/(p->W1*p->DZP[KM1]*p->DZN[KP])
                        + sigxyz2/(p->W1*p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -1.0/(p->W1*p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -1.0/(p->W1*p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -1.0/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -1.0/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -sigxyz2/(p->W1*p->DZP[KM1]*p->DZN[KP])     
                        - p->sigxx[FIJK]/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -sigxyz2/(p->W1*p->DZP[KM1]*p->DZN[KM1]) 
                        + p->sigxx[FIJK]/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
            
            
            d->rhsvec.V[n] +=  2.0*p->sigx[FIJK]*(PSI[FIp1JKp1] - PSI[FIm1JKp1] - PSI[FIp1JKm1] + PSI[FIm1JKm1])
                            /(p->W1*(p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + 2.0*p->sigy[FIJK]*(PSI[FIJp1Kp1] - PSI[FIJm1Kp1] - PSI[FIJp1Km1] + PSI[FIJm1Km1])
                            /(p->W1*(p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        }
        
        if(p->wet[IJ]==0 || p->deep[IJ]==0 || p->flag7[FIJK]<0)
        {
        d->M.p[n]  =  1.0;


        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

        d->M.w[n] = 0.0;
        d->M.e[n] = 0.0;

        d->M.t[n] = 0.0;
        d->M.b[n] = 0.0;
        
        d->rhsvec.V[n] =  0.0;
        }
	
	++n;
	}
    
    
    
    
    
    /*
    
    
    
    
    n=0;
	LOOP
	{
        if((p->X10==0 || a->fb(i,j,k)>0.0) && p->flagsf4[IJK]>0 && (a->phi(i,j,k)>=0.0 || p->I21==0))
        {
            
		if((p->flag4[Im1JK]<0 && bc(i-1,j,k)==0) || (p->X10==1 && a->fb(i-1,j,k)<0.0)
           || (p->flagsf4[Im1JK]<0 && bc(i-1,j,k)==0) || (a->phi(i-1,j,k)<0.0 && p->I21==1))
		{
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
        
        if(p->flag4[Im1JK]<0 && bc(i-1,j,k)==1)
		{
        a->rhsvec.V[n] += a->M.s[n]*p->Ui*p->DXP[IM1];
		a->M.p[n] += a->M.s[n];
		a->M.s[n] = 0.0;
		}
		
		if((p->flag4[Ip1JK]<0 && bc(i+1,j,k)==0) || (p->X10==1 && a->fb(i+1,j,k)<0.0)
           || (p->flagsf4[Ip1JK]<0 && bc(i+1,j,k)==0) || (a->phi(i+1,j,k)<0.0 && p->I21==1))
		{
		a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
        
        if(p->flag4[Ip1JK]<0 && bc(i+1,j,k)==2)
		{
         a->rhsvec.V[n] -= a->M.n[n]*p->Uo*p->DXP[IP1];
         a->M.p[n] += a->M.n[n];
		a->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0 || (p->X10==1 && a->fb(i,j-1,k)<0.0)
           || p->flagsf4[IJm1K]<0 || (a->phi(i,j-1,k)<0.0 && p->I21==1))
		{
		a->M.p[n] += a->M.e[n];
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0 || (p->X10==1 && a->fb(i,j+1,k)<0.0)
           || p->flagsf4[IJp1K]<0 || (a->phi(i,j+1,k)<0.0 && p->I21==1))
		{
		a->M.p[n] += a->M.w[n];
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0 || (p->X10==1 && a->fb(i,j,k-1)<0.0)
           || p->flagsf4[IJKm1]<0 || (a->phi(i,j,k-1)<0.0 && p->I21==1))
		{
		a->M.p[n] += a->M.b[n];
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0 || (p->X10==1 && a->fb(i,j,k+1)<0.0)
           || p->flagsf4[IJKp1]<0 || (a->phi(i,j,k+1)<0.0 && p->I21==1))
		{
		a->M.p[n] += a->M.t[n];
		a->M.t[n] = 0.0;
		}
        }

	++n;
	}*/
    
}
