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

#include"ioflow_void.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"fdm_nhf.h"
#include"vrans_v.h"
#include"vrans_f.h"
#include"rheology_v.h"
#include"rheology_f.h"
#include"turbulence.h"
#include"patchBC_interface.h"

ioflow_v::ioflow_v(lexer *p, ghostcell *pgc, patchBC_interface *ppBC)  : flowfile_in(p,pgc)
{
    pBC = ppBC;
    
	tanphi=0.0;
    if(p->W101>0)
    tanphi=tan(p->W102_phi*(PI/180.0));
}

ioflow_v::~ioflow_v()
{
}

void ioflow_v::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void ioflow_v::discharge(lexer *p, fdm* a, ghostcell* pgc)
{
    // patchBC
    pBC->patchBC_discharge(p,a,pgc);
}

void ioflow_v::inflow(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    if(p->I230>0)
    ff_inflow(p,a,pgc,u,v,w);
    
    prheo->filltau(p,a,pgc);
    
    velocity_inlet(p,a,pgc,u,v,w);
    
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_v::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    velocity_inlet(p,a,pgc,u,v,w);
    
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_v::velocity_inlet(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
//  velocity inlet
    GC1LOOP
    {
        if(p->W11==1)
        if(p->gcb1[n][3]==1 && p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        u(i-1,j,k) = p->W11_u;
        u(i-2,j,k) = p->W11_u;
        u(i-3,j,k) = p->W11_u;
        }
        
        if(p->W12==1)
        if(p->gcb1[n][3]==2 && p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        u(i,j+1,k) = p->W12_u;
        u(i,j+2,k) = p->W12_u;
        u(i,j+3,k) = p->W12_u;
        }
        
        if(p->W13==3)
        if(p->gcb1[n][3]==3 && p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        u(i,j-1,k) = p->W13_u;
        u(i,j-2,k) = p->W13_u;
        u(i,j-3,k) = p->W13_u;
        }
        
        if(p->W14==1)
        if(p->gcb1[n][3]==4 && p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        u(i+1,j,k) = p->W14_u;
        u(i+2,j,k) = p->W14_u;
        u(i+3,j,k) = p->W14_u;
        }
        
        if(p->W15==1)
        if(p->gcb1[n][3]==5 && p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        u(i,j,k-1) = p->W15_u;
        u(i,j,k-2) = p->W15_u;
        u(i,j,k-3) = p->W15_u;
        }
        
        if(p->W16==1)
        if(p->gcb1[n][3]==6 && p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        u(i,j,k) = p->W16_u;
        u(i,j,k+1) = p->W16_u;
        u(i,j,k+2) = p->W16_u;
        u(i,j,k+3) = p->W16_u;
        }
    }
    
    GC2LOOP
    {
        if(p->W11==1)
        if(p->gcb2[n][3]==1 && p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        v(i-1,j,k) = p->W11_v;
        v(i-2,j,k) = p->W11_v;
        v(i-3,j,k) = p->W11_v;
        }
        
        if(p->W12==1)
        if(p->gcb2[n][3]==2 && p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        v(i,j+1,k) = p->W12_v;
        v(i,j+2,k) = p->W12_v;
        v(i,j+3,k) = p->W12_v;
        }
        
        if(p->W13==3)
        if(p->gcb2[n][3]==3 && p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        v(i,j-1,k) = p->W13_v;
        v(i,j-2,k) = p->W13_v;
        v(i,j-3,k) = p->W13_v;
        }
        
        if(p->W14==1)
        if(p->gcb2[n][3]==4 && p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        v(i+1,j,k) = p->W14_v;
        v(i+2,j,k) = p->W14_v;
        v(i+3,j,k) = p->W14_v;
        }
        
        if(p->W15==1)
        if(p->gcb2[n][3]==5 && p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        v(i,j,k-1) = p->W15_v;
        v(i,j,k-2) = p->W15_v;
        v(i,j,k-3) = p->W15_v;
        }
        
        if(p->W16==1)
        if(p->gcb2[n][3]==6 && p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        v(i,j,k+1) = p->W16_v;
        v(i,j,k+2) = p->W16_v;
        v(i,j,k+3) = p->W16_v;
        }
        
    }
    
    
    
    
    
    GC3LOOP
    {
        if(p->W11==1)
        if(p->gcb3[n][3]==1 && p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        w(i-1,j,k) = p->W11_w;
        w(i-2,j,k) = p->W11_w;
        w(i-3,j,k) = p->W11_w;
        }
        
        if(p->W12==1)
        if(p->gcb3[n][3]==2 && p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        w(i,j+1,k) = p->W12_w;
        w(i,j+2,k) = p->W12_w;
        w(i,j+3,k) = p->W12_w;
        }
        
        if(p->W13==3)
        if(p->gcb3[n][3]==3 && p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        w(i,j-1,k) = p->W13_w;
        w(i,j-2,k) = p->W13_w;
        w(i,j-3,k) = p->W13_w;
        }
        
        if(p->W14==1)
        if(p->gcb3[n][3]==4 && p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        w(i+1,j,k) = p->W14_w;
        w(i+2,j,k) = p->W14_w;
        w(i+3,j,k) = p->W14_w;
        }
        
        if(p->W15==1)
        if(p->gcb3[n][3]==5 && p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        w(i,j,k-1) = p->W15_w;
        w(i,j,k-2) = p->W15_w;
        w(i,j,k-3) = p->W15_w;
        }
        
        if(p->W16==1)
        if(p->gcb3[n][3]==6 && p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        w(i,j,k+1) = p->W16_w;
        w(i,j,k+2) = p->W16_w;
        w(i,j,k+3) = p->W16_w;
        }
        
    }
}

void ioflow_v::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->I230>0)
    ff_waterlevel(p,a,pgc,a->phi);
    
    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void ioflow_v::fsfrkout(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
}

void ioflow_v::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
    pBC->patchBC_waterlevel(p,a,pgc,f);
}

void ioflow_v::fsfrkoutV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_v::fsfrkinV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_v::fsfrkoutVa(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_v::fsfrkinVa(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_v::iogcb_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void  ioflow_v::isource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
    double porousterm;

	count=0;
    if(p->B240>0 && p->B241==1)
    ULOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->u(i,j,k) + 0.5*p->B240_C[n]*a->u(i,j,k)*fabs(a->u(i,j,k));
		}
	
    a->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
	pvrans->u_source(p,a);
    
    //Rheology
    prheo->u_source(p,a);
}

void  ioflow_v::jsource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
    double porousterm;

	count=0;
    if(p->B240>0 && p->B242==1)
    VLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->v(i,j,k) + 0.5*p->B240_C[n]*a->v(i,j,k)*fabs(a->v(i,j,k));
		}
	
    a->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
	pvrans->v_source(p,a);
    
    //Rheology
    prheo->v_source(p,a);
}

void  ioflow_v::ksource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
    double porousterm;
	
	count=0;
    if(p->B240>0 && p->B243==1)
    WLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*a->visc(i,j,k)*a->w(i,j,k) + 0.5*p->B240_C[n]*a->w(i,j,k)*fabs(a->w(i,j,k));
		}

    a->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
	pvrans->w_source(p,a);
    
    //Rheology
    prheo->w_source(p,a);
}

void ioflow_v::isource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	d->rhsvec.V[n]=0.0;
	
    double porousterm;

	// Darcy Porosity
	count=0;
    if(p->B240>0 && p->B241==1)
    LOOP
	{
		
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*d->VISC[IJK]*d->U[IJK] + 0.5*p->B240_C[n]*d->U[IJK]*fabs(d->U[IJK]);
		}
	
    d->rhsvec.V[count] -= porousterm;
	++count;
	}
	
	//VRANS
   //pvrans->u_source(p,a);
}

void ioflow_v::jsource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	d->rhsvec.V[n]=0.0;
	
    double porousterm;

	count=0;
    if(p->B240>0 && p->B242==1)
    VLOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*d->VISC[IJK]*d->V[IJK] + 0.5*p->B240_C[n]*d->V[IJK]*fabs(d->V[IJK]);
		}
	
    d->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
    //pvrans->v_source(p,a);
}

void ioflow_v::ksource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	d->rhsvec.V[n]=0.0;
	
    double porousterm;
	
	count=0;
    if(p->B240>0 && p->B243==1)
    LOOP
	{
		// porous media
		porousterm=0.0;
		for(n=0;n<p->B240;++n)
		{
			if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
			if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
			if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
			porousterm=p->B240_D[n]*d->VISC[IJK]*d->W[IJK] + 0.5*p->B240_C[n]*d->W[IJK]*fabs(d->W[IJK]);
		}

    d->rhsvec.V[count] -= porousterm;
	++count;
	}
    
    //VRANS
    //pvrans->w_source(p,a);
}

void ioflow_v::pressure_io(lexer *p, fdm *a, ghostcell* pgc)
{
    double pval=0.0;

        GC4LOOP
        if(p->gcb4[n][4]==2)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
		pval=0.0;
		
			if(p->B77==1)
			{
			pval=(p->phiout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
			
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
		
			if(p->B77==2)
			{
			double eps,H;
                
            eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
            if(a->phi(i,j,k)>eps)
            H=1.0;

            if(a->phi(i,j,k)<-eps)
            H=0.0;

            if(fabs(a->phi(i,j,k))<=eps)
            H=0.5*(1.0 + a->phi(i,j,k)/eps + (1.0/PI)*sin((PI*a->phi(i,j,k))/eps));
        
            pval=(1.0-H)*a->press(i,j,k);
            
             a->press(i,j,k)=pval;
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
			
        }
        
    pBC->patchBC_pressure(p,a,pgc,a->press);
}

void ioflow_v::turbulence_io(lexer *p, fdm* a, ghostcell* pgc)
{
}

void ioflow_v::u_relax(lexer *p, fdm *a, ghostcell *pgc, field &uvel)
{
	double epsi,H,fbval;
    double dist;
    double cosb,sinb;

    for(int qn=0; qn<p->W41; ++qn)
    {
        cosb = cos(p->W41_beta[qn]*PI/180.0);
        
        ULOOP
        {
            dist = sqrt(pow(p->W41_xc[qn]-p->pos1_x(),2.0) + pow(p->W41_yc[qn]-p->pos_y(),2.0));
            
            if(dist>epsi)
            H=1.0;

            if(dist<-epsi)
            H=0.0;

            if(fabs(dist)<=epsi)
            H=0.5*(1.0 + dist/epsi + (1.0/PI)*sin((PI*dist)/epsi));	
            
        
            
            if(0.5*(a->phi(i,j,k)+a->phi(i+1,j,k))>0.0)
            a->u(i,j,k) = H*a->u(i,j,k) + (1.0-H)*p->W41_vel[qn]*cosb;
        }
    }
    
}

void ioflow_v::v_relax(lexer *p, fdm *a, ghostcell *pgc, field &vvel)
{
	double epsi,H,fbval;
    double dist;
    double cosb,sinb;
    
    for(int qn=0; qn<p->W41; ++qn)
    {
        sinb = sin(p->W41_beta[qn]*PI/180.0);
        
        VLOOP
        {
            dist = sqrt(pow(p->W41_xc[qn]-p->pos_x(),2.0) + pow(p->W41_yc[qn]-p->pos2_y(),2.0));
            
            if(dist>epsi)
            H=1.0;

            if(dist<-epsi)
            H=0.0;

            if(fabs(dist)<=epsi)
            H=0.5*(1.0 + dist/epsi + (1.0/PI)*sin((PI*dist)/epsi));	
            
        
            
            if(0.5*(a->phi(i,j,k)+a->phi(i,j+1,k))>0.0)
            a->v(i,j,k) = H*a->v(i,j,k) + (1.0-H)*p->W41_vel[qn]*sinb;
        }
    }
}

void ioflow_v::w_relax(lexer *p, fdm *a, ghostcell *pgc, field &wvel)
{
	
}

void ioflow_v::p_relax(lexer *p, fdm *a, ghostcell *pgc, field &press)
{
    /*double tau0,tau,pval,phival,H,gamma;
    double epsi = 1.6*p->DXM;
    
    if(p->W1
    LOOP
    {
        phival = a->phi(i,j,k);
    
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));   
    
    
        // get gamma from rheology
        
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0=tanphi*pval + p->W102_c;
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*pval + p->W102_c)*(1.0-exp(-p->W103*gamma));
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_u also use new input W 104

        if(p->count==0)
        tau0=p->W96;
        
        // get tau from rheology
        
        if(tau<tau0)
        a->press(i,j,k) = H*a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22) + (1.0-H)*a->press(i,j,k);
        
        
    }*/
}

void ioflow_v::phi_relax(lexer *p, ghostcell *pgc, field &f)
{
}

void ioflow_v::vof_relax(lexer *p, ghostcell *pgc, field &f)
{
}

void ioflow_v::turb_relax(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
}

void ioflow_v::U_relax(lexer *p, ghostcell *pgc, double *U, double *UH)
{
}

void ioflow_v::V_relax(lexer *p, ghostcell *pgc, double *V, double *VH)
{
}

void ioflow_v::W_relax(lexer *p, ghostcell *pgc, double *W, double *WH)
{
}

void ioflow_v::P_relax(lexer *p, ghostcell *pgc, double *P)
{
}

void ioflow_v::WL_relax(lexer *p, ghostcell *pgc, slice &WL, slice &depth)
{
}

void ioflow_v::fi_relax(lexer *p, ghostcell *pgc, field &f, field &phi)
{
}

void ioflow_v::fivec_relax(lexer *p, ghostcell *pgc, double *f)
{
}

void ioflow_v::fifsf_relax(lexer *p, ghostcell *pgc, slice& f)
{
}

void ioflow_v::visc_relax(lexer *p, ghostcell *pgc, slice& f)
{
}

void ioflow_v::eta_relax(lexer *p, ghostcell *pgc, slice &f)
{
}

void ioflow_v::um_relax(lexer *p, ghostcell *pgc, slice &P, slice &bed, slice &eta)
{
}

void ioflow_v::vm_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
}

void ioflow_v::wm_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
}

void ioflow_v::ws_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
}

void ioflow_v::pm_relax(lexer *p, ghostcell *pgc, slice &f)
{
}

double ioflow_v::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    double val=0.0;

    return val;
}

int ioflow_v::iozonecheck(lexer *p, fdm*a)
{	
	int check =1;
	
	return check;
}

void ioflow_v::inflow_walldist(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, reini *preini, ioflow *pflow)
{
}

void ioflow_v::discharge2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // patchBC
    pBC->patchBC_discharge2D(p,b,pgc,b->P,b->Q,b->eta,b->bed);
}

void ioflow_v::Qin2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
}

void ioflow_v::Qout2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
}

void ioflow_v::inflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_v::rkinflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_v::isource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
	SLICELOOP1
	b->F(i,j)=0.0;
}

void ioflow_v::jsource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
	SLICELOOP2
	b->G(i,j)=0.0;
}

void ioflow_v::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    if(p->B269==0)
	pvrans = new vrans_v(p,pgc);
	
	if(p->B269==1 || p->S10==2)
	pvrans = new vrans_f(p,pgc);
    
    if(p->W90==0)
    prheo = new rheology_v(p,a);
    
    if(p->W90==1)
    prheo = new rheology_f(p,a);
}

void ioflow_v::full_initialize2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void ioflow_v::flowfile(lexer *p, fdm* a, ghostcell* pgc, turbulence *pturb)
{
}

void ioflow_v::wavegen_precalc(lexer *p, ghostcell *pgc)
{
    
}

void ioflow_v::wavegen_precalc_ini(lexer *p, ghostcell *pgc)
{
    
}

void ioflow_v::wavegen_2D_precalc(lexer *p, fdm2D *b, ghostcell *pgc)
{
    
}

void ioflow_v::wavegen_2D_precalc_ini(lexer *p, ghostcell *pgc)
{
    
}

void ioflow_v::ini2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void ioflow_v::ini_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void ioflow_v::ini_ptf(lexer *p, fdm* a, ghostcell* pgc)
{
    
}

void ioflow_v::veltimesave(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans) 
{
    
}

void ioflow_v::inflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, double *Fi, double *Uin,slice &Fifsf, slice &eta)
{

}

void ioflow_v::rkinflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, slice &frk, slice &f)
{
}

void ioflow_v::vrans_sed_update(lexer *p,fdm *a,ghostcell *pgc,vrans *pvrans)
{
    pvrans->sed_update(p,a,pgc);
}

void ioflow_v::ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{

}

void ioflow_v::discharge_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{

}

void ioflow_v::inflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{

}

void ioflow_v::rkinflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{

}

void ioflow_v::wavegen_precalc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}

void ioflow_v::wavegen_precalc_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}

void ioflow_v::waterlevel2D(lexer *p, fdm2D *b, ghostcell* pgc, slice &eta)
{
    pBC->patchBC_waterlevel2D(p,b,pgc,eta);
}

void ioflow_v::fsfinflow_nhflow(lexer *p, fdm_nhf* d, ghostcell* pgc, slice &WL)
{

}