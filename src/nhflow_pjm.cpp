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

#include"nhflow_pjm.h"
#include"lexer.h"
#include"nhflow_curv.h"
#include"fdm_nhf.h"
#include"nhflow_membrane_beta.h"
#include"nhflow_thinbody.h"
#include"nhflow_bodyflux.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"nhflow_pjm_ho.h"
#include"solver.h"
#include"ioflow.h"
#include"density_f.h"
#include"patchBC_interface.h"
#include"vrans.h"
#include"vrans_definitions.h"

#define WLVL (fabs(WL(i,j))>(p->A544)?WL(i,j):1.0e20)

nhflow_pjm::nhflow_pjm(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(1.0)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson(p);

    gcval_press=540;
    
    gamma=0.5;
    
    p->Darray(P0,p->imax*p->jmax*(p->kmax+2));
}

nhflow_pjm::~nhflow_pjm()
{
    delete pd;
    delete ppois;
    delete [] P0;
}

void nhflow_pjm::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL,
                        double *UH, double *VH, double *WH, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
    
    FLOOP
    P0[FIJK] = d->P[FIJK];

    rhs(p,d,pgc,d->U,d->V,d->W,alpha);

    ppois->start(p,d,d->P);

        starttime=pgc->timer();

    psolv->startF(p,pgc,d->P,d->rhsvec,d->M,8);

        endtime=pgc->timer();

	pgc->start7P(p,d->P,gcval_press);
    
    // cut faces (X 15 10): the pressure of the nodes without water from their neighbours
    if(d->thinbody!=nullptr && d->thinbody->cutmode)
    d->thinbody->extend_pressure(p,d,pgc,d->P);

	ucorr(p,d,WL,UH,d->P,alpha);
	vcorr(p,d,WL,VH,d->P,alpha);
	wcorr(p,d,WL,WH,d->P,alpha);
    
    // membranes (X 330): further projection passes. The collocated correction (wide gradient) and the
    // compact matrix differ most where the mobility jumps; repeated passes converge the wide divergence
    // of the corrected velocity to zero (the compact matrix acts as preconditioner).
    for(int it=1; it<d->MPROJ; ++it)
    {
        nhflow_membrane_velupdate(p,d,pgc,WL,UH,VH,WH);
        
        rhs(p,d,pgc,d->U,d->V,d->W,alpha);
        
        FLOOP
        P0[FIJK] = 0.0;
        
        ppois->start(p,d,P0);
        psolv->startF(p,pgc,P0,d->rhsvec,d->M,8);
        pgc->start7P(p,P0,gcval_press);
        
        // the corrections read d->P: apply the increment, then add it to the pressure
        double *Psave = d->P;
        d->P = P0;
        ucorr(p,d,WL,UH,d->P,alpha);
        vcorr(p,d,WL,VH,d->P,alpha);
        wcorr(p,d,WL,WH,d->P,alpha);
        d->P = Psave;
        
        FLOOP
        d->P[FIJK] += P0[FIJK];
        
        pgc->start7P(p,d->P,gcval_press);
    }
    
    // membranes (X 330): face corrections for the Rhie-Chow continuity flux of the next stage
    if(d->MBETA!=nullptr)
    nhflow_membrane_rc_store(p,d,pgc,d->P,alpha*p->dt);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

// mesh refinement: start() without the solve (the membrane passes X 330 are not supported there; the sharp hull,
// X 15 8, makes one pass)
double* nhflow_pjm::amr_prepare(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha)
{
    FLOOP
    P0[FIJK] = d->P[FIJK];
    
    rhs(p,d,pgc,d->U,d->V,d->W,alpha);
    ppois->start(p,d,d->P);
    
    return d->P;
}

void nhflow_pjm::amr_finish(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *UH, double *VH, double *WH, double alpha)
{
	pgc->start7P(p,d->P,gcval_press);
    
	ucorr(p,d,WL,UH,d->P,alpha);
	vcorr(p,d,WL,VH,d->P,alpha);
	wcorr(p,d,WL,WH,d->P,alpha);
    
    // the hull as a sharp thin body (X 15 8, one projection pass): face corrections for the Rhie-Chow continuity
    // flux of the next stage, as start()
    if(d->MBETA!=nullptr)
    nhflow_membrane_rc_store(p,d,pgc,d->P,alpha*p->dt);
}

void nhflow_pjm::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *UH, double *P, double alpha)
{
    const bool ho = nhflow_pjm_ho::active(p,d);
    
    // cut faces (X 15 10): the horizontal pressure force of the cells as finite volumes (nhflow_cut_gx)
    const bool cut = d->thinbody!=nullptr && d->thinbody->cutmode;

    // membranes (X 330): face mobilities in the horizontal gradient, see nhflow_membrane_beta.h
    // cut faces (X 15 10): the correction of a cut cell per its water (nhflow_thinbody::ik)
    // curvilinear grid: the horizontal gradient by the chain rule (nhflow_curv.h)
	LOOP
    WETDRYDEEP
    {
    double px=0.0, py=0.0;
    if(p->geo_curv==1)
    nhflow_curv::pgrad_layer(p,d->P,i,j,k,px,py);
    
	UH[IJK] -= alpha*p->dt*CPORNH*WL(i,j)*(1.0/p->W1)*(d->thinbody!=nullptr ? d->thinbody->ik(IJK) : 1.0)*
    
                (cut ? nhflow_cut_gx(p,d,WL,d->P,i,j,k) :
                (d->MBETA!=nullptr ? nhflow_membrane_gradx(p,d,d->P,i,j,k) :
                p->geo_curv==1 ? px :
                (ho && nhflow_pjm_ho::x(p,d,i,j,k)) ? nhflow_pjm_ho::gradx(p,d->P,i,j,k) :     // fourth order (nhflow_pjm_ho.h)
                ((0.5*(d->P[FIp1JKp1]+d->P[FIp1JK])-0.5*(d->P[FIm1JKp1]+d->P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1])))
                
                + MBETAVAL*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *VH, double *P, double alpha)
{
    const bool ho = nhflow_pjm_ho::active(p,d);
    
    // cut faces (X 15 10): the horizontal pressure force of the cells as finite volumes (nhflow_cut_gx)
    const bool cut = d->thinbody!=nullptr && d->thinbody->cutmode;

    if(p->j_dir==1)
    LOOP
    WETDRYDEEP
    {
    double px=0.0, py=0.0;
    if(p->geo_curv==1)
    nhflow_curv::pgrad_layer(p,d->P,i,j,k,px,py);
    
    VH[IJK] -= alpha*p->dt*CPORNH*WL(i,j)*(1.0/p->W1)*(d->thinbody!=nullptr ? d->thinbody->ik(IJK) : 1.0)*
    
                (cut ? nhflow_cut_gy(p,d,WL,d->P,i,j,k) :
                (d->MBETA!=nullptr ? nhflow_membrane_grady(p,d,d->P,i,j,k) :
                p->geo_curv==1 ? py :
                (ho && nhflow_pjm_ho::y(p,d,i,j,k)) ? nhflow_pjm_ho::grady(p,d->P,i,j,k) :
                ((0.5*(d->P[FIJp1Kp1]+d->P[FIJp1K])-0.5*(d->P[FIJm1Kp1]+d->P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1])))
                
                + MBETAVAL*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *WH, double *P, double alpha)
{
    LOOP
    WETDRYDEEP
	WH[IJK] -= MBETAVAL*alpha*p->dt*CPORNH*(1.0/p->W1)*(d->thinbody!=nullptr ? d->thinbody->ik(IJK) : 1.0)*((d->P[FIJKp1]-d->P[FIJK])/(p->DZN[KP]));
}

void nhflow_pjm::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    double U1,U2,V1,V2,fac;
    double dUdz,dVdz;
    double dWdz;
    
    n=0;
    FBASELOOP
    {
	d->rhsvec.V[n]=0.0;
    ++n;
    }

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
    
    // sharp thin bodies (X 330 'mobility sharp'): wall velocity at the blocked links
    if(d->thinbody!=nullptr)
    d->thinbody->projection_rhs(p,d,U,V,W,alpha);
    
    // direct forcing with the body in the projection (X 15 11): the flux of the body through the closed part of its faces
    if(nhflow_bodyflux::active(p,d))
    nhflow_bodyflux::rhs(p,d,U,V,W,d->FBUX,d->FBUY,d->FBWZ,alpha*p->dt,true);
}

void nhflow_pjm::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    // curvilinear grid: g eta sum(d_face S_face)/A (nhflow_curv.h)
    if(p->geo_curv==1)
    LOOP
    WETDRY
    d->F[IJK] += nhflow_curv::bedslope(p,d,1,i,j);
    
    // X 13 2: the level of the hydrostatic terms (d->zeta), eta otherwise
    slice4 &E = d->zeta!=nullptr ? *d->zeta : d->eta;
    
    if(p->geo_curv==0)
    LOOP
    WETDRY
    d->F[IJK] += E(i,j)*fabs(p->W22)*
                (d->dfx(i,j) - d->dfx(i-1,j))/(p->DXN[IP]);
}

void nhflow_pjm::vpgrad(lexer*p,fdm_nhf *d, slice &WL)
{
    if(p->geo_curv==1)
    LOOP
    WETDRY
    d->G[IJK] += nhflow_curv::bedslope(p,d,2,i,j);
    
    slice4 &E = d->zeta!=nullptr ? *d->zeta : d->eta;
    
    if(p->geo_curv==0)
    LOOP
    WETDRY
	d->G[IJK] += E(i,j)*fabs(p->W22)*
                 (d->dfy(i,j) - d->dfy(i,j-1))/(p->DYN[JP]);
}

void nhflow_pjm::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
}

