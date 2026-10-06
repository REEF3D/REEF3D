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
    p->Darray(WLN,p->imax*p->jmax);
    wl_ini=0;
}

nhflow_suspended_IM1::~nhflow_suspended_IM1()
{
}

void nhflow_suspended_IM1::start(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_scalar_convection *pconvec, nhflow_diffusion *pdiff, solver *psolv, ioflow *pflow, sediment_fdm *s)
{
    starttime=pgc->timer();
    drysave(p,d,s);
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
    
    // depth belonging to the new concentration, for the D^n/D^(n+1) ratio of the next time step
    SLICELOOP4
    WLN[IJ] = d->WL(i,j);
    wl_ini=1;
    
	p->susptime=pgc->timer()-starttime;
	p->suspiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"suspiter: "<<p->suspiter<<"  susptime: "<<setprecision(3)<<p->susptime<<endl;
}

void nhflow_suspended_IM1::timesource(lexer* p, fdm_nhf *d, double *FN)
{
    // conservative form for D*C with the transport coefficients divided by D^(n+1):
    //   (D^(n+1) C^(n+1) - D^n C^n)/(D^(n+1) dt) = C^(n+1)/dt - (D^n/D^(n+1)) C^n/dt
    // D^n: depth at the end of the last solve (includes bed changes and the flow step since then)
    int count=0;
    double hn,ho;

    LOOP
    {
        hn = MAX(d->WL(i,j),1.0e-20);
        ho = wl_ini==1?MAX(WLN[IJ],0.0):hn;
        
        d->M.p[count]+= 1.0/p->dt;

        d->rhsvec.V[count] += d->L[IJK] + (ho/hn)*d->CONC[IJK]/p->dt;

	++count;
    }
}

void nhflow_suspended_IM1::ctimesave(lexer *p, fdm_nhf *d)
{
}

void nhflow_suspended_IM1::drysave(lexer *p, fdm_nhf *d, sediment_fdm *s)
{
    // columns that fell dry since the last solve: the solve sets their concentration to zero,
    // so the sediment they still hold (D^n C^n) is handed to the bed (deposited by the Exner step)
    if(wl_ini==0)
    return;
    
    SLICELOOP4
    if(p->wet[IJ]==0)
    {
        KLOOP
        PCHECK
        s->dryd(i,j) += MAX(d->CONC[IJK],0.0)*p->DZN[KP]*MAX(WLN[IJ],0.0);
    }
}

void nhflow_suspended_IM1::fill_wvel(lexer *p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s)
{
    // WVEL: vertical volume flux across the sigma faces for the conservative (form=2) ifou scheme:
    // omegaF is the face volume flux D*dsigma/dt [m/s],
    // settling across a sigma face is -ws (D*dsigma/dt of -ws = -ws).
    // Face k lies between cells k-1 and k; bed (k=0) and surface (k=knoz) faces stay closed,
    // bed exchange is handled by suspsource().
    double ws_eff,nval,Re_p,cface;

    Re_p = s->ws*p->S20/p->W2;
    nval = (4.7 + 0.41*pow(Re_p,0.75))/(1 + 0.175*pow(Re_p,0.75));

    FLOOP
    {
    WVEL[FIJK] = 0.0;

        if(k>0 && k<p->knoz && p->DF[IJK]>0 && p->wet[IJ]==1)
        {
        cface = 0.5*(d->CONC[IJK] + d->CONC[IJKm1]);
        ws_eff = s->ws * pow(MAX(1.0 - cface/0.635, 0.0), nval);
        WVEL[FIJK] = d->omegaF[FIJK] - ws_eff;
        }
    }
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
        d->rhsvec.V[count]  += (-s->ws)*(-s->cbe(i,j))/zdist;
        d->M.p[count] += (s->ws)/zdist;
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
            // closed faces (walls, domain edges, solid and dry neighbours, bed, free surface):
            // implicit zero gradient, the off-diagonal is folded into the diagonal, so the
            // diffusive flux through the face is exactly zero (explicit C^n left a flux ~ C^(n+1)-C^n);
            // at inflow edges the inflow concentration is the cell value
            if(p->flag4[Im1JK]<0 || p->DF[Im1JK]<0 || p->wet[Im1J]==0)
            {
            d->M.p[n] += d->M.s[n];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0 || p->wet[Ip1J]==0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0 || p->wet[IJm1]==0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0 || p->wet[IJp1]==0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            
            // bed: no diffusive flux, the exchange with the bed is in suspsource (implicit zero gradient)
            if(k==0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKm1]<0 || p->DF[IJKm1]<0) && k>0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k<p->knoz-1)
            {
            d->M.p[n] += d->M.t[n];
            d->M.t[n] = 0.0;
            }
            
            // free surface: impermeable, implicit zero gradient (was a ghost value C = 0: diffusive loss)
            if(k==p->knoz-1)
            {
            d->M.p[n] += d->M.t[n];
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
            s->cb(i,j) = MAX(MIN(d->CONC[IJK],0.1),0.0);

            if(p->S61==2)
            s->cb(i,j) = Rouse_formula(p,d,s,d->CONC[IJK]);
        }
    }    
    pgc->gcsl_start4(p,s->cb,1);
    
    
    double Uh;
    
    SEDSLICELOOP
    s->qbs(i,j) = 0.0;
    
    SEDSLICELOOP
    KLOOP
    {
        Uh = sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK]);
        
        s->qbs(i,j) += Uh*d->CONC[IJK]*p->DZN[KP]*d->WL(i,j);
    }
    
    pgc->gcsl_start4(p,s->qbs,1);
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
