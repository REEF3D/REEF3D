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
    p->Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    
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
    
    ini_bc(p,d,pgc);

    starttime=pgc->timer();
	
	pgc->start49V(p,PSI,gcval_pot);
	
	LOOP
	PSI[IJK] = 0.0;
    
    pgc->start49V(p,PSI,gcval_pot);
	
    itermem=p->N46;
    p->N46=5000;
	
    laplace(p,d,pgc);
    psolv->startV(p,pgc,PSI,d->rhsvec,d->M,44);
    
    pgc->start49V(p,PSI,gcval_pot);
    
    ucalc(p,d);
	vcalc(p,d);
	wcalc(p,d);
    
	pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
    
    pgc->start4V(p,d->UH,14);
    pgc->start4V(p,d->VH,15);
    pgc->start4V(p,d->WH,16);
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
    
    pgc->start4V(p,d->UH,14);
    pgc->start4V(p,d->VH,15);
    pgc->start4V(p,d->WH,16);

    endtime=pgc->timer();
    p->laplaceiter=p->solveriter;
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0  && (p->count%p->P12==0))
	cout<<"lapltime: "<<p->laplacetime<<"  lapiter: "<<p->laplaceiter<<endl<<endl;

    p->N46=itermem;
    

    finalize:

    p->del_Iarray(BC,p->imax*p->jmax*(p->kmax+3));
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
    
    
    LOOP
    if(p->DF[IJK]<0 || p->DF[Im1JK]<0 || p->DF[Ip1JK]<0)
    {
	d->U[IJK] = 0.0;
    d->UH[IJK] = 0.0;
    }
}

void nhflow_potential_f::vcalc(lexer *p, fdm_nhf *d)
{	
    LOOP
    {
	d->V[IJK] =  (PSI[IJp1K]-PSI[IJm1K])/(p->DYP[JP]+p->DYP[JP1]) 
    
              + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KM1]));
              
    d->VH[IJK] = d->WL(i,j)*d->V[IJK];
    }


    LOOP
    if(p->DF[IJK]<0 || p->DF[IJm1K]<0 || p->DF[IJp1K]<0)
    {
	d->V[IJK] = 0.0;
    d->VH[IJK] = 0.0;
    }
}

void nhflow_potential_f::wcalc(lexer *p, fdm_nhf *d)
{
    LOOP
    {
	d->W[IJK] =  p->sigz[IJ]*(PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KP1]); 
              
    d->WH[IJK] = d->WL(i,j)*d->W[IJK];
    }
    
    
	LOOP
    if(p->DF[IJK]<0 || p->DF[IJm1K]<0 || p->DF[IJp1K]<0)
    {
	d->W[IJK] = 0.0;
    d->WH[IJK] = 0.0;
    }
}

void nhflow_potential_f::rhs(lexer *p, fdm_nhf *d)
{
    count=0;
    LOOP
    {
    d->rhsvec.V[count] = 0.0;
    count++;
    }
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
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            d->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP])
                        + 1.0/(p->DXP[IM1]*p->DXN[IP])
                        
                        + 1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + 1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + sigxyz2/(p->DZP[KM1]*p->DZN[KP])
                        + sigxyz2/(p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -sigxyz2/(p->DZP[KM1]*p->DZN[KP])     
                        - p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) 
                        + p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
            
            
            d->rhsvec.V[n] =  0.0;(p->sigx[FIJK]+p->sigx[FIJKp1])*(PSI[Ip1JKp1] - PSI[Im1JKp1] - PSI[Ip1JKm1] + PSI[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + (p->sigx[FIJK]+p->sigx[FIJKp1])*(PSI[IJp1Kp1] - PSI[IJm1Kp1] - PSI[IJp1Km1] + PSI[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        }
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0)
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
    
    
  double denom,ab;  
    
    
    n=0;
	LOOP
	{
        if(p->DF[IJK]>0 && p->flag4[IJK]>0)
        {
           
		if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && BC[Im1JK]==0)
		{
		d->M.p[n] += d->M.s[n];
		d->M.s[n] = 0.0;
		}
        
        if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && BC[Im1JK]==1)
		{
        d->rhsvec.V[n] += d->M.s[n]*p->Ui*p->DXP[IM1];
		d->M.p[n] += d->M.s[n];
		d->M.s[n] = 0.0;
		}
		
		if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && BC[Ip1JK]==0)
		{
		d->M.p[n] += d->M.n[n];
		d->M.n[n] = 0.0;
		}
        
        if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && BC[Ip1JK]==2)
		{
         d->rhsvec.V[n] -= d->M.n[n]*p->Uo*p->DXP[IP1];
         d->M.p[n] += d->M.n[n];
		d->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
		{
		d->M.p[n] += d->M.e[n];
		d->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
		{
		d->M.p[n] += d->M.w[n];
		d->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]>0 && p->DF[IJKm1]<0)
		{
		d->M.p[n] += d->M.b[n];
		d->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
		{
		d->M.p[n] += d->M.t[n];
		d->M.t[n] = 0.0;
		}
        
        
         // KBEDBC
            if(p->flag4[IJKm1]<0)
            {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            ab = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + d->Bx(i,j)*p->sigx[FIJK] + d->By(i,j)*p->sigy[FIJK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    d->M.n[n] +=  ab*2.0*p->DZN[KP]*(d->Bx(i,j))/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    d->M.s[n] += -ab*2.0*p->DZN[KP]*(d->Bx(i,j))/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJp1]==1 && p->wet[IJm1]==1)
                    {
                    d->M.w[n] +=  ab*2.0*p->DZN[KP]*d->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    d->M.e[n] += -ab*2.0*p->DZN[KP]*d->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }

                d->M.t[n] += ab;
                d->M.b[n] = 0.0;
            }
        }

	++n;
	}
    
}

void nhflow_potential_f::ini_bc(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    BASELOOP
    BC[IJK]=0;
    
    LOOP
    {
        if(p->flag4[Im1JK]<0)
		BC[Im1JK]=0;
		
		if(p->flag4[Ip1JK]<0)
		BC[Ip1JK]=0;
		
		if(p->flag4[IJm1K]<0)
		BC[IJm1K]=0;
		
		if(p->flag4[IJp1K]<0)
		BC[IJp1K]=0;
		
		if(p->flag4[IJKm1]<0)
		BC[IJKm1]=0;
		
		if(p->flag4[IJKp1]<0)
		BC[IJKp1]=0;
    }
    

    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
            i=p->gcb4[n][0]; 
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];  
            
            WETDRY
            {
            if(p->gcb4[n][3]==1)
            BC[Im1JK]=1;
            
            if(p->gcb4[n][3]==3)
            BC[IJm1K]=1;
            
            if(p->gcb4[n][3]==2)
            BC[IJp1K]=1;
            
            if(p->gcb4[n][3]==4)
            BC[Ip1JK]=1;
            }
 
        }
        
        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
            i=p->gcb4[n][0]; 
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];  
            
            WETDRY
            {
            if(p->gcb4[n][3]==1)
            BC[Im1JK]=2;
            
            if(p->gcb4[n][3]==3)
            BC[IJm1K]=2;
            
            if(p->gcb4[n][3]==2)
            BC[IJp1K]=2;
            
            if(p->gcb4[n][3]==4)
            BC[Ip1JK]=2;
            }
 
        }
    }
}
