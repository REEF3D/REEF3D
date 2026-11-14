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

#include"nhflow_pjm_corr.h"
#include"lexer.h"
#include"nhflow_curv.h"
#include"fdm_nhf.h"
#include"nhflow_membrane_beta.h"
#include"ghostcell.h"
#include"nhflow_poisson_pcorr.h"
#include"nhflow_pjm_ho.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson.h"
#include"density_f.h"
#include"patchBC_interface.h"
#include"vrans.h"
#include"vrans_definitions.h"
#include"nhflow_thinbody.h"
#include"nhflow_bodyflux.h"

#define WLVL (fabs(WL(i,j))>(p->A544)?WL(i,j):1.0e20)

nhflow_pjm_corr::nhflow_pjm_corr(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson_pcorr(p);
    
    p->Darray(PCORR,p->imax*p->jmax*(p->kmax+2));
    
    gcval_press=540;
    
    gamma=0.5;
}

nhflow_pjm_corr::~nhflow_pjm_corr()
{
    delete pd;
    delete ppois;
    delete [] PCORR;
}

void nhflow_pjm_corr::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL,
                        double *UH, double *VH, double *WH, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
    
        starttime=pgc->timer();

    rhs(p,d,pgc,d->U,d->V,d->W,alpha);

    ppois->start(p,d,PCORR);
    
        p->matrixtime+=pgc->timer()-starttime;


        starttime=pgc->timer();

    psolv->startF(p,pgc,PCORR,d->rhsvec,d->M,8);

        endtime=pgc->timer();
        
    
    presscorr(p,d,pgc,WL,d->P,PCORR,alpha);
    
    pgc->start7P(p,d->P,gcval_press);
    pgc->start7P(p,PCORR,gcval_press);
    
	ucorr(p,d,WL,UH,PCORR,alpha);
	vcorr(p,d,WL,VH,PCORR,alpha);
	wcorr(p,d,WL,WH,PCORR,alpha);
    
    // membranes (X 330): further projection passes, see nhflow_pjm
    for(int it=1; it<d->MPROJ; ++it)
    {
        nhflow_membrane_velupdate(p,d,pgc,WL,UH,VH,WH);
        
        rhs(p,d,pgc,d->U,d->V,d->W,alpha);
        ppois->start(p,d,PCORR);
        psolv->startF(p,pgc,PCORR,d->rhsvec,d->M,8);
        presscorr(p,d,pgc,WL,d->P,PCORR,alpha);
        pgc->start7P(p,d->P,gcval_press);
        pgc->start7P(p,PCORR,gcval_press);
        ucorr(p,d,WL,UH,PCORR,alpha);
        vcorr(p,d,WL,VH,PCORR,alpha);
        wcorr(p,d,WL,WH,PCORR,alpha);
    }
    
    // membranes (X 330): face corrections of the pressure increment for the Rhie-Chow continuity flux.
    // The old pressure P^n acts on the cell velocities through the predictor (wide gradient), and the Poisson
    // right-hand side is the divergence of those cell velocities averaged to the faces. The face field the
    // compact matrix makes divergence free is therefore avg(U) + avg(dU(PCORR)) - dU_f(PCORR): only the
    // increment is corrected. With the total pressure the face flux contains avg(G_wide P^n) - G_compact P^n,
    // a divergence the projection never sees; it moves the free surface of the layer columns in an odd-even
    // pattern that grows (seen with moving membranes, R_t = R_n and floorpressure 0).
    // sharp thin bodies (the hull X 15 8, nhflow_thinbody): the predictor takes the old pressure with the link mobilities
    // (upgrad, vpgrad, wpgrad), so the cell velocities carry the face-mobility gradient of the total pressure and the
    // Rhie-Chow faces are stored with it
    if(d->MBETA!=nullptr)
    nhflow_membrane_rc_store(p,d,pgc,(d->thinbody!=nullptr || (p->X15==11 && p->X10>0)) ? d->P : PCORR,alpha*p->dt);

    p->poissoniter=p->solveriter;

	p->ptime=endtime-starttime;
    p->poissontime+=p->ptime;    

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->ptime<<endl;
}

// mesh refinement: start() without the solve (the membrane passes X 330 are not supported there)
double* nhflow_pjm_corr::amr_prepare(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha)
{
    rhs(p,d,pgc,d->U,d->V,d->W,alpha);
    ppois->start(p,d,PCORR);
    
    return PCORR;
}

void nhflow_pjm_corr::amr_finish(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *UH, double *VH, double *WH, double alpha)
{
    presscorr(p,d,pgc,WL,d->P,PCORR,alpha);
    
    pgc->start7P(p,d->P,gcval_press);
    pgc->start7P(p,PCORR,gcval_press);
    
	ucorr(p,d,WL,UH,PCORR,alpha);
	vcorr(p,d,WL,VH,PCORR,alpha);
	wcorr(p,d,WL,WH,PCORR,alpha);
}

void nhflow_pjm_corr::presscorr(lexer* p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *P, double *PCORR, double alpha)
{
	FLOOP
    WETDRYDEEP
    P[FIJK] += PCORR[FIJK];
}

void nhflow_pjm_corr::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    double U1,U2,V1,V2,fac;
    double dUdz,dVdz;
    double dWdz;
    
    n=0;
    FBASELOOP
    {
	d->rhsvec.V[n]=0.0;
    PCORR[FIJK]=0.0;
    ++n;
    }
    
    pgc->start7P(p,PCORR,gcval_press);
    
    const bool ho = nhflow_pjm_ho::active(p,d);
    
    n=0;
    LOOP
    {
        WETDRYDEEP
        {
        // bed node (k==0): the velocity below the bed is taken as 0 here, whatever the bed ghost cells of the velocity
        // (A 518): the bed face then carries (1-fac) U of the bed cell, the half cell of the bed node half its
        // horizontal divergence. With the mirrored ghost cells of the wall-function bed (A 518 1) in this place the
        // dispersion of the waves changed (linear sloshing, kh = pi/2, 5 layers: period +1.8 % instead of -0.1 %).
        // This is the continuity equation of the half cell on the bed, on a sloping bed too: the slope enters through
        // sigx dU/dz only. (A 522, retired, gave the bed cell a larger weight on slopes, 1 - 0.5 max(1 - 5 |slope|, 0):
        // on the Beji & Battjes bar the higher harmonics behind the bar came out wrong with it.)
        const bool bednode = (k==0);
        const double Ub1 = bednode ? 0.0 : U[Im1JKm1], Ub2 = bednode ? 0.0 : U[Ip1JKm1];
        const double Vb1 = bednode ? 0.0 : V[IJm1Km1], Vb2 = bednode ? 0.0 : V[IJp1Km1];
        const double Ubc = bednode ? 0.0 : U[IJKm1],   Vbc = bednode ? 0.0 : V[IJKm1];
        
        fac = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);    

        U1 = (1.0-fac)*U[Im1JK] + fac*Ub1; 
        U2 = (1.0-fac)*U[Ip1JK] + fac*Ub2; 
        
        V1 = (1.0-fac)*V[IJm1K] + fac*Vb1; 
        V2 = (1.0-fac)*V[IJp1K] + fac*Vb2;

        dUdz = (U[IJK] - Ubc)/p->DZP[KM1];
        dVdz = (V[IJK] - Vbc)/p->DZP[KM1];
        dWdz = p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZP[KM1];
         
        // curvilinear grid: the horizontal divergence by the chain rule (nhflow_curv.h)
        if(p->geo_curv==1)
        d->rhsvec.V[n] =      -  (nhflow_curv::divh_level(p,d,U,V,i,j,k)
                                + p->sigx[FIJK]*dUdz
                                + p->sigy[FIJK]*dVdz
                                + dWdz)/(alpha*p->dt);
        else
        {
        // fourth-order horizontal divergence (nhflow_pjm_ho.h)
        double dUdx = (U2-U1)/(p->DXP[IP] + p->DXP[IM1]);
        double dVdy = (V2-V1)/(p->DYP[JP] + p->DYP[JM1]);
        
        if(ho)
        {
        const double f0 = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
        
        if(nhflow_pjm_ho::x(p,d,i,j,k))
        dUdx = nhflow_pjm_ho::divx(p,U,i,j,k,f0);
        
        if(nhflow_pjm_ho::y(p,d,i,j,k))
        dVdy = nhflow_pjm_ho::divy(p,V,i,j,k,f0);
        }
        
        d->rhsvec.V[n] =      -  (dUdx
                                + p->sigx[FIJK]*dUdz
                                
                                + dVdy
                                + p->sigy[FIJK]*dVdz

                                + dWdz)/(alpha*p->dt);
        }
        }
                            
    ++n;
    }
    
    // sharp thin bodies: wall velocity at the blocked links (as nhflow_pjm::rhs)
    if(d->thinbody!=nullptr)
    d->thinbody->projection_rhs(p,d,U,V,W,alpha);
    
    // direct forcing with the body in the projection (X 15 11): the flux of the body through the closed part of its faces
    if(nhflow_bodyflux::active(p,d))
    nhflow_bodyflux::rhs(p,d,U,V,W,d->FBUX,d->FBUY,d->FBWZ,alpha*p->dt,true);
}

void nhflow_pjm_corr::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *UH, double *PCORR, double alpha)
{
    const bool ho = nhflow_pjm_ho::active(p,d);
    
	LOOP
    WETDRYDEEP
    {
    // membranes (X 330): face mobilities in the horizontal gradient, see nhflow_membrane_beta.h
    if(d->MBETA!=nullptr)
    dPdx = nhflow_membrane_gradx(p,d,PCORR,i,j,k);
    else
    if(p->geo_curv==1)
    nhflow_curv::pgrad_layer(p,PCORR,i,j,k,dPdx,dPdy);
    else
    dPdx = (0.5*(PCORR[FIp1JKp1]+PCORR[FIp1JK])-0.5*(PCORR[FIm1JKp1]+PCORR[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]);
    
    // fourth-order horizontal gradient (nhflow_pjm_ho.h)
    if(ho && nhflow_pjm_ho::x(p,d,i,j,k))
    dPdx = nhflow_pjm_ho::gradx(p,PCORR,i,j,k);
    
	UH[IJK] -= alpha*p->dt*CPORNH*WL(i,j)*(1.0/p->W1)*
    
                (dPdx
                
                + MBETAVAL*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((PCORR[FIJKp1]-PCORR[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_corr::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *VH, double *PCORR, double alpha)
{
    const bool ho = nhflow_pjm_ho::active(p,d);
    
    if(p->j_dir==1)
    LOOP
    WETDRYDEEP
    {
    if(d->MBETA!=nullptr)
    dPdy = nhflow_membrane_grady(p,d,PCORR,i,j,k);
    else
    if(p->geo_curv==1)
    nhflow_curv::pgrad_layer(p,PCORR,i,j,k,dPdx,dPdy);
    else
    dPdy = (0.5*(PCORR[FIJp1Kp1]+PCORR[FIJp1K])-0.5*(PCORR[FIJm1Kp1]+PCORR[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]);
    
    if(ho && nhflow_pjm_ho::y(p,d,i,j,k))
    dPdy = nhflow_pjm_ho::grady(p,PCORR,i,j,k);
    
    VH[IJK] -= alpha*p->dt*CPORNH*WL(i,j)*(1.0/p->W1)*
    
                (dPdy
                
                + MBETAVAL*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((PCORR[FIJKp1]-PCORR[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_corr::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *WH, double *PCORR, double alpha)
{
    LOOP
    WETDRYDEEP
	WH[IJK] -= MBETAVAL*alpha*p->dt*CPORNH*(1.0/p->W1)*((PCORR[FIJKp1]-PCORR[FIJK])/(p->DZN[KP]));
}

void nhflow_pjm_corr::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    const bool ho = nhflow_pjm_ho::active(p,d);
    
    // curvilinear grid: g eta sum(d_face S_face)/A (nhflow_curv.h)
    if(p->geo_curv==1)
    LOOP
    WETDRY
    d->F[IJK] += nhflow_curv::bedslope(p,d,1,i,j);
    
    if(p->geo_curv==0)
    LOOP
    WETDRY
    d->F[IJK] += (d->zeta!=nullptr ? (*d->zeta)(i,j) : d->eta(i,j))*fabs(p->W22)*
                (d->dfx(i,j) - d->dfx(i-1,j))/(p->DXN[IP]);
               
    LOOP
    WETDRYDEEP
    {
    if(p->geo_curv==1)
    nhflow_curv::pgrad_layer(p,d->P,i,j,k,dPdx,dPdy);
    else
    dPdx = (0.5*(d->P[FIp1JKp1]+d->P[FIp1JK])-0.5*(d->P[FIm1JKp1]+d->P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]);
    
    // sharp thin bodies: the old pressure with the link mobilities, as the correction (no pressure across the hull)
    const bool tb = (d->thinbody!=nullptr || (p->X15==11 && p->X10>0)) && d->MBETA!=nullptr;   // X 15 11: the open fractions of the links
    
    if(tb)
    dPdx = nhflow_membrane_gradx(p,d,d->P,i,j,k);
    
    // fourth-order horizontal gradient (nhflow_pjm_ho.h)
    if(ho && nhflow_pjm_ho::x(p,d,i,j,k))
    dPdx = nhflow_pjm_ho::gradx(p,d->P,i,j,k);
    
    d->F[IJK] -= (1.0/p->W1)*WL(i,j)*
                (dPdx
                + (tb ? MBETAVAL : 1.0)*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_corr::vpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    const bool ho = nhflow_pjm_ho::active(p,d);
    
    if(p->geo_curv==1)
    LOOP
    WETDRY
    d->G[IJK] += nhflow_curv::bedslope(p,d,2,i,j);
    
    if(p->j_dir==1 && p->geo_curv==0)
    LOOP
    WETDRY
	d->G[IJK] += (d->zeta!=nullptr ? (*d->zeta)(i,j) : d->eta(i,j))*fabs(p->W22)*
                 (d->dfy(i,j) - d->dfy(i,j-1))/(p->DYN[JP]);
       
    if(p->j_dir==1) 
    LOOP
    WETDRYDEEP
    {
    if(p->geo_curv==1)
    nhflow_curv::pgrad_layer(p,d->P,i,j,k,dPdx,dPdy);
    else
    dPdy = (0.5*(d->P[FIJp1Kp1]+d->P[FIJp1K])-0.5*(d->P[FIJm1Kp1]+d->P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]);
    
    const bool tb = (d->thinbody!=nullptr || (p->X15==11 && p->X10>0)) && d->MBETA!=nullptr;   // X 15 11: the open fractions of the links
    
    if(tb)
    dPdy = nhflow_membrane_grady(p,d,d->P,i,j,k);
    
    if(ho && nhflow_pjm_ho::y(p,d,i,j,k))
    dPdy = nhflow_pjm_ho::grady(p,d->P,i,j,k);
    
	d->G[IJK] -= (1.0/p->W1)*WL(i,j)*
                (dPdy
                + (tb ? MBETAVAL : 1.0)*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_corr::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRYDEEP
    d->H[IJK] -= ((d->thinbody!=nullptr || (p->X15==11 && p->X10>0)) && d->MBETA!=nullptr ? MBETAVAL : 1.0)*(1.0/p->W1)*((d->P[FIJKp1]-d->P[FIJK])/(p->DZN[KP]));
}

