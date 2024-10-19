/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"suspended_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"sediment_fdm.h"

suspended_RK3::suspended_RK3(lexer* p, fdm* a) : wvel(p)
{
	gcval_susp=60;
}

suspended_RK3::~suspended_RK3()
{
}

void suspended_RK3::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow, sediment_fdm *s)
{
    field4 ark1(p),ark2(p);
    fill_wvel(p,a,pgc,s);    bcsusp_start(p,a,pgc,s,a->conc);
    
// Step 1
    starttime=pgc->timer();    clearrhs(p,a);
    suspsource(p,a,a->conc,s);
    pconvec->start(p,a,a->conc,4,a->u,a->v,wvel);
	pdiff->diff_scalar(p,a,pgc,psolv,a->conc,a->visc,a->eddyv,1.0,1.0);

	LOOP
	ark1(i,j,k) = a->conc(i,j,k)
                + p->dt*a->L(i,j,k);
	
    bcsusp_start(p,a,pgc,s,ark1);
    sedfsf(p,a,ark1);
	pgc->start4(p,ark1,gcval_susp);

// Step 2    clearrhs(p,a);
    suspsource(p,a,ark1,s);
    pconvec->start(p,a,ark1,4,a->u,a->v,wvel);
	pdiff->diff_scalar(p,a,pgc,psolv,ark1,a->visc,a->eddyv,1.0,0.25);

	LOOP
	ark2(i,j,k) = 0.75*a->conc(i,j,k)
                + 0.25*ark1(i,j,k)
                + 0.25*p->dt*a->L(i,j,k);
	
    bcsusp_start(p,a,pgc,s,ark2);
    sedfsf(p,a,ark2);
	pgc->start4(p,ark2,gcval_susp);

// Step 3    clearrhs(p,a);
    suspsource(p,a,ark2,s);
    pconvec->start(p,a,ark2,4,a->u,a->v,wvel);
	pdiff->diff_scalar(p,a,pgc,psolv,ark2,a->visc,a->eddyv,1.0,2.0/3.0);

	LOOP
	a->conc(i,j,k) = (1.0/3.0)*a->conc(i,j,k)
				  + (2.0/3.0)*ark2(i,j,k)
				  + (2.0/3.0)*p->dt*a->L(i,j,k);
	
    bcsusp_start(p,a,pgc,s,a->conc);
    sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);    fillconc(p,a,s);

	p->susptime=pgc->timer()-starttime;
}

void suspended_RK3::ctimesave(lexer *p, fdm* a)
{
}

void suspended_RK3::fill_wvel(lexer *p, fdm* a, ghostcell *pgc, sediment_fdm *s)
{
    WLOOP
    wvel(i,j,k) = a->w(i,j,k) - s->ws;
    
    pgc->start3(p,wvel,12);
}

void suspended_RK3::suspsource(lexer* p,fdm* a,field& conc, sediment_fdm *s)
{
    LOOP
    a->L(i,j,k)=0.0;
    
    LOOP
    if(p->flagsf4[IJK]<0)
    conc(i,j,k) = 0.0;
        
    GCDF4LOOP
    {
            i=p->gcdf4[n][0];
            j=p->gcdf4[n][1];
            k=p->gcdf4[n][2];
    
        if(s->cbe(i,j)>=conc(i,j,k))
        a->L(i,j,k) += s->ws/p->DZN[KP]*(s->cbe(i,j)-conc(i,j,k));
        
       //if(s->cbe(i,j)<conc(i,j,k))
        //a->L(i,j,k) += s->ws/p->DZN[KP]*s->cbe(i,j);
        
        s->cb(i,j)=conc(i,j,k);
    }

}

void suspended_RK3::bcsusp_start(lexer* p, fdm* a,ghostcell *pgc, sediment_fdm *s, field& conc)
{

        LOOP
        if(p->flagsf4[IJK]<0)
        conc(i,j,k) = 0.0;
            
        GCDF4LOOP
        {
            i=p->gcdf4[n][0];
            j=p->gcdf4[n][1];
            k=p->gcdf4[n][2];
            
            //conc(i,j,k) =    s->cb(i,j);
            conc(i,j,k-1) =  conc(i,j,k);
            conc(i,j,k-2) =  conc(i,j,k);
            conc(i,j,k-3) =  conc(i,j,k);
        }
    
    
    // Inflow
    /*double cval;
    double zdist =1.0;
    double zcoor;
    double adist;
    
    for(n=0;n<p->gcin_count;++n)
    if(p->gcin[n][3]>0)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    adist=0.05*s->waterlevel(i,j);
    
    zcoor = a->topo(i,j,k) + 0.000001;
    
    cval=0.0;
    
    if(a->topo(i,j,k)>0.0 && a->phi(i,j,k)>0.0)
    {
    cval = 0.001*pow(((s->waterlevel(i,j)-zcoor)/zcoor)*(adist/(s->waterlevel(i,j)-adist + 0.000001)),1.0);
    
    //cout<<"cval: "<<setprecision(6)<<cval<<" zwl: "<<setprecision(6)<<s->waterlevel(i,j)<<" zcoor: "<<setprecision(6)<<zcoor<<endl;
    }
    
    conc(i-1,j,k) = cval;
    conc(i-2,j,k) = cval;
    conc(i-3,j,k) = cval;
    }*/
    }
void suspended_RK3::fillconc(lexer* p, fdm* a, sediment_fdm *s){    double dist;    double d50=p->S20;    double adist=0.5*d50;    double deltab=3.0*d50;    double cx,cy;        if(p->S34==1)    GC4LOOP    if(p->gcb4[n][4]==5)    {        i=p->gcb4[n][0];        j=p->gcb4[n][1];        k=p->gcb4[n][2];                s->conc(i,j) = a->conc(i,j,k);                //dist = 0.5*p->DZN[KP]-adist;                //s->conc(i,j) = (s->cbe(i,j)*(dist-deltab+adist) + a->conc(i,j,k+1)*(deltab-adist))/(dist);        //if(s->conc(i,j)>s->cbe(i,j))        //cout<<"conc: "<<s->conc(i,j)<<" cbe: "<<s->cbe(i,j)<<endl;    }
    
    if(p->S34==1)
    GCDF4LOOP
    {
        i=p->gcdf4[n][0];
        j=p->gcdf4[n][1];
        k=p->gcdf4[n][2];
        
        s->conc(i,j) = a->conc(i,j,k);
        
        //dist = p->ZP[KP1]-s->bedzh(i,j)-adist;
        
        //s->conc(i,j) = (s->cbe(i,j)*(dist-deltab+adist) + a->conc(i,j,k+1)*(deltab-adist))/(dist);
        
        //if(s->conc(i,j)>s->cbe(i,j))
        //cout<<"conc: "<<s->conc(i,j)<<" cbe: "<<s->cbe(i,j)<<endl;
    }        if(p->S34==2)    ILOOP    JLOOP    {        cx=0.0;        cy=0.0;            KLOOP        PCHECK        {        cx += 0.5*(a->u(i,j,k) + a->u(i-1,j,k))*a->conc(i,j,k)*p->DZN[KP];        cy += 0.5*(a->v(i,j,k) + a->v(i,j-1,k))*a->conc(i,j,k)*p->DZN[KP];        }    s->conc(i,j) = sqrt(cx*cx + cy*cy);    }}
void suspended_RK3::sedfsf(lexer* p,fdm* a,field& conc)
{
    LOOP
    if(a->phi(i,j,k)<0.0)
    conc(i,j,k)=0.0;
}

void suspended_RK3::clearrhs(lexer* p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	++count;
    }
}