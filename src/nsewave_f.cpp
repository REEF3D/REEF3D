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

#include"nsewave_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"
#include"momentum.h"
#include"sflow_eta_weno.h"
#include"sflow_hxy_weno.h"

nsewave_f::nsewave_f(lexer *p, fdm *a, ghostcell *pgc, heat *&pheat, concentration *&pconc) : 
                depth(p),bed(p),L(p),hp(p),hx(p),hy(p)
{
	// bed ini
	SLICELOOP4
	bed(i,j) = p->bed[IJ];
	pgc->gcsl_start4(p,bed,50);
	
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

    pupdate = new fluid_update_fsf(p,a,pgc);
    
    p->phimean=p->F60;
    
    LOOP
    a->phi(i,j,k) = a->eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    
    
    if(p->j_dir==0)        
    epsi = p->A440*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = p->A440*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
}

nsewave_f::~nsewave_f()
{
}

void nsewave_f::start(lexer* p, fdm* a, ghostcell* pgc, momentum *pmom, diffusion *pdiff, turbulence *pturb,
                      convection* pconvec, pressure *ppress, poisson *ppois, solver *ppoissonsolv, solver *psolv, 
                      ioflow* pflow, vrans* pvrans, sixdof *p6dof, vector<net*>& pnet)
{
    
    // Momentum
    pmom->start(p,a,pgc,pvrans,p6dof,pnet);
    
    
    
    // fill eta_n
    SLICELOOP4
	{
    a->eta_n(i,j) = a->eta(i,j);
	L(i,j)=0.0;
	}
    pgc->gcsl_start4(p,a->eta_n,gcval_phi);
    
    
    // Calculate Eta
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

	
    IULOOP
	JULOOP
    {
		KULOOP
        UFLUIDCHECK
		{
		phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
		
			if(phival>epsi)
			H=1.0;

			if(phival<-epsi)
			H=0.0;

			if(fabs(phival)<=epsi)
			H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
			
			a->P(i,j) += a->u(i,j,k)*p->DZN[KP]*H;
		}
    }
    
    IVLOOP
	JVLOOP
	{	
			KVLOOP
            VFLUIDCHECK
			{
			 phival = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));
			
				if(phival>epsi)
				H=1.0;

				if(phival<-epsi)
				H=0.0;

				if(fabs(phival)<=epsi)
				H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
				
				a->Q(i,j) += a->v(i,j,k)*p->DZN[KP]*H;

			}
    }
	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
	

    SLICELOOP4
    a->eta(i,j) -= p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);	  

    
    pflow->eta_relax(p,pgc,a->eta);
    
    pgc->gcsl_start4(p,a->eta,gcval_phi);
    
    FLUIDLOOP
    a->phi(i,j,k) = a->eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    
    pupdate->start(p,a,pgc);
    
}

void nsewave_f::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
    pupdate->start(p,a,pgc);
}

void nsewave_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
    p->phimean=p->F60;
    
    pflow->eta_relax(p,pgc,a->eta);
    
    pgc->gcsl_start4(p,a->eta,gcval_phi);
    
    FLUIDLOOP
    a->phi(i,j,k) = a->eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);
    
    pupdate->start(p,a,pgc);
}


