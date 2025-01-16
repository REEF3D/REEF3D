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
Author: Tobias Martin, Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"initialize.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"
#include"interpolation.h"

VOF_PLIC::VOF_PLIC
(
    lexer* p,
    fdm *a,
    ghostcell* pgc,
    heat *pheat
):gradient(p),norm_vec(p),alpha(p),nx(p),ny(p),nz(p),vof1(p),vof2(p),vof3(p),phival(p),V_w_p(p),V_w_m(p),alphastore(p),phistep(p),phiS0(p),phiS1(p),phiS2(p),vofstep(p),vofS0(p),vofS1(p),vofS2(p),
    phiaux(p),Vx_p(p),Vx_m(p),Vz_p(p),Vz_m(p),Vn_p(p),Vn_m(p),F_x(p),F_z(p),F_n(p),F_new(p),Flux_x(p),Flux_z(p),Crossflux_xz(p),Crossflux_zx(p)
{
    if(p->F50==1)
    gcval_frac=71;

    if(p->F50==2)
    gcval_frac=72;

    if(p->F50==3)
    gcval_frac=73;

    if(p->F50==4)
    gcval_frac=74;

    pupdate = new fluid_update_vof(p,a,pgc);

    reini_ = new reini_RK3(p,1);
    
    ipol = new interpolation(p);
    
    if(p->j_dir>0)
        Sweepdim=3;
    else
        Sweepdim=2;

    sSweep = -1;
    
    ininorVecLS(p);
    
    S_S[0][0]=0;
    S_S[0][1]=1;
    S_S[0][2]=2;
    S_S[1][0]=2;
    S_S[1][1]=1;
    S_S[1][2]=0;
    S_S[2][0]=1;
    S_S[2][1]=0;
    S_S[2][2]=2;
    S_S[3][0]=0;
    S_S[3][1]=2;
    S_S[3][2]=1;
    S_S[4][0]=2;
    S_S[4][1]=0;
    S_S[4][2]=1;
    S_S[5][0]=1;
    S_S[5][1]=2;
    S_S[5][1]=0;
    
    S_2D[0][0]=0;
    S_2D[0][1]=2;
    S_2D[1][0]=2;
    S_2D[1][1]=0;
}

VOF_PLIC::~VOF_PLIC()
{
}


void VOF_PLIC::update
(

    lexer *p,
    fdm *a,
    ghostcell *pgc,
    field &F
)
{
    pupdate->start(p,a,pgc);
}

void VOF_PLIC::start
(
    fdm* a,
    lexer* p,
    convection* pconvec,
    solver* psolv,
    ghostcell* pgc,
    ioflow* pflow,
    reini* preini,
    particle_corr* ppls,
    field &F
)
{
    
    if(p->count<=10 && p->F89==1)
    {
        reini_->start(a,p,a->phi,pgc,pflow);
        reini_->start(a,p,a->phi,pgc,pflow);
        reini_->start(a,p,a->phi,pgc,pflow);
    }
   // pflow->fsfinflow(p,a,pgc);

   // pflow->vof_relax(p,a,pgc,a->vof);
    pgc->start4(p,a->phi,1);
    pgc->start4(p,a->vof,1);
    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);
    pgc->start3(p,a->w,12);
    starttime=pgc->timer();
    
    if(p->j_dir>0)
    {
        if(sSweep<5)
            sSweep++;
        else
            sSweep=0;
    }
    else
    {
        if(sSweep<1)
            sSweep++;
        else
            sSweep=0;
    }
    if(p->F90<=3)
        stepwise_scheme(a,p,pgc);
    else if(p->F90==4);
        symmetric_scheme2D(a,p,pgc);
    pgc->start4(p,vofstep,1);
    LOOP
    {   
        
        
        if(vofstep(i,j,k)<0.0 ) //|| ((p->F89>0 && p->F91>0.0) && phistep(i,j,k)< -p->F91*p->psi))
            vofstep(i,j,k)=0.0;
        if(vofstep(i,j,k)>1.0 ) //|| ((p->F89>0 && p->F91>0.0) && phistep(i,j,k)>p->F91*p->psi))
            vofstep(i,j,k)=1.0;
            
     //   sprayfilter(a,p);
            
        if(vofstep(i,j,k)>0.001 && vofstep(i,j,k)<0.999)
        {   
            reconstructPlane_alt(a,p,vofstep);
        }
        else
        {    
            nx(i,j,k)=2.0;
            ny(i,j,k)=0.0;
            nz(i,j,k)=2.0;
            alpha(i,j,k)=1E20;
        }
        
        phiaux(i,j,k)=1E05;
        a->vof(i,j,k)=vofstep(i,j,k);
    }
    pgc->start4(p,nx,1);
    pgc->start4(p,ny,1);
    pgc->start4(p,nz,1);
    pgc->start4(p,alpha,1);
    pgc->start4(p,phiaux,1);
    pgc->start4(p,vofstep,1);
    pgc->start4(p,a->vof,1);
        LOOP
        {
            if((vofstep(i,j,k)>0.001 && vofstep(i,j,k)<0.999))
            {
                redistancePhiByPlane_Bonn(a,p);
            }
        }
    
        pgc->start4(p,phiaux,1);
        LOOP
        {
            if(fabs(phiaux(i,j,k))<1E04)
            {
                a->phi(i,j,k)=phiaux(i,j,k);
            }
            else if(p->F89==1)
            {
                a->phi(i,j,k)=phistep(i,j,k);
            }
            else
            {
                if(a->vof(i,j,k)>=0.5)
                    a->phi(i,j,k)=0.5;
                else
                    a->phi(i,j,k)=-0.5;
            }
               
        }
    
        pgc->start4(p,a->phi,1);
  
        if(p->F89==1)
            reini_->start(a,p,a->phi,pgc,pflow);
    
        pgc->start4(p,a->phi,1);

    pupdate->start(p,a,pgc);
   
   
}

