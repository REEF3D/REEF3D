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

#include"nhflow_suspended_IM1.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_scalar_convection.h"
#include"nhflow_diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"sediment_fdm.h"

nhflow_suspended_IM1::nhflow_suspended_IM1(lexer* p) 
{
	gcval_susp=60;

    p->Darray(WVEL,p->imax*p->jmax*(p->kmax+2));
}

nhflow_suspended_IM1::~nhflow_suspended_IM1()
{
}

void nhflow_suspended_IM1::start(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_scalar_convection *pconvec, nhflow_diffusion *pdiff, solver *psolv, ioflow *pflow, sediment_fdm *s)
{
    starttime=pgc->timer();
    clearrhs(p,d);
    fill_wvel(p,d,pgc,s);
    pconvec->start(p,d,d->CONC,4,d->U,d->V,WVEL);
    pdiff->diff_scalar(p,d,pgc,psolv,d->CONC,1.0,1.0);
	suspsource(p,d,d->CONC,s);
	timesource(p,d,d->CONC);
    bcsusp_start(p,d,pgc,s,d->CONC);
    psolv->startV(p,pgc,d->CONC,d->rhsvec,d->M,4);
	pgc->start60V(p,d->CONC,gcval_susp);
    fillconc(p,d,pgc,s);
	p->susptime=pgc->timer()-starttime;
	p->suspiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"suspiter: "<<p->suspiter<<"  susptime: "<<setprecision(3)<<p->susptime<<endl;
}

void nhflow_suspended_IM1::timesource(lexer* p, fdm_nhf *d, double *FN)
{
    int count=0;

    LOOP
    {
        d->M.p[count]+= 1.0/p->dt;

        d->rhsvec.V[count] += d->L[IJK] + d->CONC[IJK]/p->dt;

	++count;
    }
}

void nhflow_suspended_IM1::ctimesave(lexer *p, fdm_nhf *d)
{
}

void nhflow_suspended_IM1::fill_wvel(lexer *p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s)
{
    LOOP
    {
    WVEL[IJK] = 0.0;
    
    if(p->DF[IJK]>0 && p->wet[IJ]==1)
    WVEL[IJK] = d->W[IJK] - s->ws;
    }
    
    pgc->start4V(p,WVEL,12);
}

void nhflow_suspended_IM1::suspsource(lexer* p, fdm_nhf *d, double *CONC, sediment_fdm *s)
{    
    double zdist;
    
    count=0;
    LOOP
    {   
        if(k==0 && p->DF[IJK]>0 && p->wet[IJ]==1)
        {
        zdist = p->DZN[KP]*d->WL(i,j);
        d->rhsvec.V[count]  += (-s->ws)*(s->cb(i,j)-s->cbe(i,j))/zdist;
        }
        
        /*
        if(p->mpirank==0)
        if(i==10 && k==p->knoz-1)
        d->rhsvec.V[count] += 0.00001;*/

	++count;
    }
}

void nhflow_suspended_IM1::bcsusp_start(lexer *p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s, double *CONC)
{
    double cval;
    
        n=0;
        LOOP
        {
            if(p->DF[IJK]>0 && p->wet[IJ]==1)
            {
                
            if(p->flag4[Im1JK]<0 || p->DF[Im1JK]<0 || p->wet[Im1J]==0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*CONC[IJK];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0 || p->wet[Ip1J]==0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*CONC[IJK];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0 || p->wet[IJm1]==0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*CONC[IJK];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0 || p->wet[IJp1]==0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*CONC[IJK];
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*CONC[IJK];
            d->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k<p->knoz-1)
            {
            d->rhsvec.V[n] -= d->M.t[n]*CONC[IJK];
            d->M.t[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k==p->knoz-1)
            {
            d->rhsvec.V[n] -= d->M.t[n]*0.0;
            d->M.t[n] = 0.0;
            }
            }

        ++n;
        }
        
        
    // turn off inside direct forcing body
        n=0;
        LOOP
        {
            if(p->DF[IJK]<0 || p->wet[IJ]==0)
            {
            d->M.p[n] = 1.0;

            d->M.n[n] = 0.0;
            d->M.s[n] = 0.0;

            d->M.w[n] = 0.0;
            d->M.e[n] = 0.0;

            d->M.t[n] = 0.0;
            d->M.b[n] = 0.0;
            
            d->rhsvec.V[n] = 0.0;
            }
        ++n;
        }
}

void nhflow_suspended_IM1::fillconc(lexer* p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s)
{
    k=0;
    SLICELOOP4
    {
    if(p->DF[IJK]<0 || p->wet[IJ]==0)
    s->cb(i,j) = 0.0;
    
        if(p->DF[IJK]>0 && p->wet[IJ]==1)
        {
        if(p->S61==1)
        s->cb(i,j) = d->CONC[IJK];
        
        if(p->S61==2)
        s->cb(i,j) = Rouse_formula(p,d,s,d->CONC[IJK]);
        }
    }    
    
    pgc->gcsl_start4(p,s->cb,1);
}

double nhflow_suspended_IM1::Rouse_formula(lexer* p, fdm_nhf *d, sediment_fdm *s, double Cc)
{
    double Ca;    
    double za,zc,P;
    
    za = 2.0*p->S20;
    
    zc = 0.5*p->DZN[KP]*p->WL[IJ];
    
    P = s->ws/(0.4* (s->shearvel_eff(i,j)>0.0?s->shearvel_eff(i,j):1.0e-6) );
    
    P = MAX(P,0.8);
    P = MIN(P,2.5);
    
    
    Ca = Cc * pow( ((p->WL[IJ]-za)/za) / ((p->WL[IJ]-zc)/zc), P);
    
    Ca = MIN(Ca,0.1);
    
    //cout<<"Cc: "<<Cc<<" Ca: "<<Ca<<" | P: "<<P<<" "<<s->shearvel_eff(i,j)<<endl;

    return Ca;
}

void nhflow_suspended_IM1::clearrhs(lexer* p, fdm_nhf *d)
{
    n=0;
    LOOP
    {    
    d->rhsvec.V[n]=0.0;
    d->L[IJK]=0.0;
    
    
            d->M.p[n] = 0.0;

            d->M.n[n] = 0.0;
            d->M.s[n] = 0.0;

            d->M.w[n] = 0.0;
            d->M.e[n] = 0.0;

            d->M.t[n] = 0.0;
            d->M.b[n] = 0.0;
	++n;
    }
}
