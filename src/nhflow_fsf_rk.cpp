/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_void.h"
#include"heat.h"
#include"concentration.h"
#include"momentum.h"
#include"sflow_hxy_weno.h"
#include"sflow_hxy_cds.h"
#include"sflow_hxy_fou.h"
#include"patchBC_interface.h"
#include"nhflow_flux_HLL.h"
#include"nhflow_flux_HLLC.h"
#include"nhflow_flux_FOU.h"

nhflow_fsf_rk::nhflow_fsf_rk(lexer *p, fdm_nhf* d, ghostcell *pgc, ioflow *pflow, patchBC_interface *ppBC) : epsi(p->A440*p->DXM),P(p),Q(p),K(p)
{
    pBC = ppBC;
    
	pupdate = new fluid_update_void();
    
    if(p->A541==1)
	phxy = new sflow_hxy_fou(p,pBC);
	
	if(p->A541==2)
	phxy = new sflow_hxy_cds(p,pBC);
	
	if(p->A541==4)
	phxy = new sflow_hxy_weno(p,pBC);
    
    if(p->A542==1)
    pfluxfsf = new nhflow_flux_FOU(p,pBC);
    
    if(p->A542==2)
    pfluxfsf = new nhflow_flux_HLL(p,pBC);
    
    if(p->A542==3)
    pfluxfsf = new nhflow_flux_HLLC(p,pBC);
    
    p->Darray(Fx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fy,p->imax*p->jmax*(p->kmax+2));
    
    wd_criterion=0.00005;
    
    if(p->A244==1)
    wd_criterion=p->A244_val;
    
    if(p->A245==1)
    wd_criterion=p->A245_val*p->DXM;
}

nhflow_fsf_rk::~nhflow_fsf_rk()
{
}

void nhflow_fsf_rk::start(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow)
{
}

void nhflow_fsf_rk::update(lexer *p, fdm_nhf* d, ghostcell *pgc, slice &f)
{
    //pupdate->start(p,a,pgc);
}

void nhflow_fsf_rk::step1(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{
    wetdry(p,d,pgc,U,V,W,d->eta);
    
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += 0.5*(U[IJK] + U[Ip1JK]);

    VLOOP
	Q(i,j) += 0.5*(V[IJK] + V[IJp1K]);
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,d->eta,P,Q);
    
    pgc->gcsl_start1(p,d->hx,10);
    pgc->gcsl_start2(p,d->hy,11);
    
    pfluxfsf->face_flux_3D(p,pgc,d,d->eta,U,V,Fx,Fy);
    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    K(i,j) += -p->DZN[KP]*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP]  + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    pgc->start4V(p,d->test,10);
    
    SLICELOOP4
    etark1(i,j) = d->eta(i,j) + p->dt*K(i,j);
     
    pflow->eta_relax(p,pgc,etark1);
    pgc->gcsl_start4(p,etark1,1);
    
    SLICELOOP4
    d->WL(i,j) = (etark1(i,j) + p->wd - d->bed(i,j));
    
    //wetdry(p,d,pgc,U,V,W,etark1);
    
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    breaking(p,d,pgc,etark1,d->eta,1.0);
}

void nhflow_fsf_rk::step2(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{
    wetdry(p,d,pgc,U,V,W,etark1);
    
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += 0.5*(U[IJK] + U[Ip1JK]);

    VLOOP
	Q(i,j) += 0.5*(V[IJK] + V[IJp1K]);
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,etark1,P,Q);
    
    pgc->gcsl_start1(p,d->hx,10);
    pgc->gcsl_start2(p,d->hy,11);
    
    pfluxfsf->face_flux_3D(p,pgc,d,etark1,U,V,Fx,Fy);
    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    K(i,j) += -p->DZN[KP]*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP]  + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    SLICELOOP4
    etark2(i,j) = 0.75*d->eta(i,j) + 0.25*etark1(i,j) + 0.25*p->dt*K(i,j);

    pflow->eta_relax(p,pgc,etark2);
    pgc->gcsl_start4(p,etark2,1);
    
    SLICELOOP4
    d->WL(i,j) =(etark2(i,j) + p->wd - d->bed(i,j));
    
    //wetdry(p,d,pgc,U,V,W,etark2);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    breaking(p,d,pgc,etark2,etark1,0.25);
}

void nhflow_fsf_rk::step3(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{
    // fill eta_n
    SLICELOOP4
    {
    d->eta_n(i,j) = d->eta(i,j);
    d->WL_n0(i,j) = d->WL_n1(i,j);
    }    
    pgc->gcsl_start4(p,d->eta_n,1);
    
    // ---
    wetdry(p,d,pgc,U,V,W,etark2);
        
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += 0.5*(U[IJK] + U[Ip1JK]);

    VLOOP
	Q(i,j) += 0.5*(V[IJK] + V[IJp1K]);
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,etark2,P,Q);
    
    pgc->gcsl_start1(p,d->hx,10);
    pgc->gcsl_start2(p,d->hy,11);
    
    pfluxfsf->face_flux_3D(p,pgc,d,etark2,U,V,Fx,Fy);
    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    K(i,j) += - p->DZN[KP]*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP]  + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    SLICELOOP4
    d->eta(i,j) = (1.0/3.0)*d->eta(i,j) + (2.0/3.0)*etark2(i,j) + (2.0/3.0)*p->dt*K(i,j);


    pflow->eta_relax(p,pgc,d->eta);
    pgc->gcsl_start4(p,d->eta,1);
    
    SLICELOOP4
    d->WL(i,j) = (d->eta(i,j) + p->wd - d->bed(i,j));
    
    SLICELOOP4
    d->WL_n1(i,j) = d->WL(i,j);
    
    //SLICELOOP4
    //d->detadt(i,j) = (d->eta(i,j) -d->eta_n(i,j))/p->dt;
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    //wetdry(p,d,pgc,U,V,W,d->eta);
    
    LOOP
    d->test[IJK] = Fx[IJK];
    
    breaking(p,d,pgc,d->eta,etark2,2.0/3.0);
}


