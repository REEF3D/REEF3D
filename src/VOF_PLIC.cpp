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
):gradient(p),norm_vec(p),alpha(p),nx(p),ny(p),nz(p),vof1(p),vof2(p),vof3(p),vof_old(p),V_w_update(p),V_a_update(p),phival(p),Watersafe(p),V_w_p_star(p),V_w_m_star(p),vof_prevstep(p),V_w_old(p),V_a_old(p),V_w_p(p),V_w_m(p),FX_p(p),FX_m(p),FZ_p(p),FZ_m(p),alphastore(p),nxstore(p),nystore(p),nzstore(p),phistep(p),phiS0(p),phiS1(p),phiS2(p),vofstep(p),vofS0(p),vofS1(p),vofS2(p),phiaux(p)
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
    
    if(p->count<=10)
    {
        reini_->start(a,p,a->phi,pgc,pflow);
        reini_->start(a,p,a->phi,pgc,pflow);
        reini_->start(a,p,a->phi,pgc,pflow);
    }
    pflow->fsfinflow(p,a,pgc);
    
    pflow->vof_relax(p,a,pgc,a->vof);
    pgc->start4(p,a->phi,1);
    pgc->start4(p,a->vof,1);
    

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
    
    int sweep;
    int Sweepdim;
    
    if(p->j_dir>0)
        Sweepdim=3;
    else
        Sweepdim=2;
        
    for(int nSweep=0 ; nSweep<Sweepdim; nSweep++)
    {
        LOOP
        {
            V_w_m(i,j,k)=0.0;
            V_w_p(i,j,k)=0.0;
        }
        pgc->start4(p,V_w_m,1);
        pgc->start4(p,V_w_p,1);
        
        if(p->j_dir>0)
            sweep=S_S[sSweep][nSweep];
        else
            sweep=S_2D[sSweep][nSweep];
        
        if(nSweep==0)
        {
            pgc->start4(p,a->phi,1);
            pgc->start4(p,a->vof,1);
            LOOP
            {
                phistep(i,j,k)=a->phi(i,j,k);
                vofstep(i,j,k)=a->vof(i,j,k);
            }
            pgc->start4(p,phistep,1);
            pgc->start4(p,vofstep,1);
        }
        else if(nSweep==1)
        {
            pgc->start4(p,phiS0,1);
            pgc->start4(p,vofS0,1);
            LOOP
            {
                phistep(i,j,k)=phiS0(i,j,k);
                vofstep(i,j,k)=vofS0(i,j,k);
            }
            pgc->start4(p,phistep,1);
            pgc->start4(p,vofstep,1);
        }
        else
        {
            pgc->start4(p,phiS1,1);
            pgc->start4(p,vofS1,1);
            LOOP
            {
                phistep(i,j,k)=phiS1(i,j,k);
                vofstep(i,j,k)=vofS1(i,j,k);
            }
            pgc->start4(p,phistep,1);
            pgc->start4(p,vofstep,1);
        }
        
        pgc->start1(p,a->u,10);
        pgc->start2(p,a->v,11);
        pgc->start3(p,a->w,12);
        
        LOOP
        {
            switch(p->F89)
                {
                    case 0:
                        break;
                    case 1:
                        transportPhi_Bonn(a,p,nSweep,sweep);
                        break;
                }
            
            bool bordercheck=false;
            
         /*   if(p->j_dir>0)
            {
                for(int isearch=-1; isearch<2; isearch++)
                {
                    for(int jsearch=-1; jsearch<2; jsearch++)
                    {
                        for(int ksearch=-1; ksearch<2; ksearch++)
                        {
                            if(phistep(i,j,k)*phistep(i+isearch,j+jsearch,k+ksearch)<0.0)
                                bordercheck=true;
                        }
                    }
                }
            }
            
            else
            {
                for(int isearch=-1; isearch<2; isearch++)
                {
                    for(int ksearch=-1; ksearch<2; ksearch++)
                    {
                        if(phistep(i,j,k)*phistep(i+isearch,j,k+ksearch)<0.0)
                                bordercheck=true;
                    }
                }
            }*/
            if((vofstep(i,j,k)>0.0001 && vofstep(i,j,k)<0.9999))
            {
                reconstructPlane_alt(a,p,vofstep(i,j,k));
                
                switch(p->F90)
                {
                    case 0:
                        break;
                    case 1:
                        advectPlane_forBonnScheme(a,p,sweep);
                        break;
                    case 3:
                        advectPlane_NewWang(a,p,sweep);
                        break;
                }
                
            }
            else if( vofstep(i,j,k)>0.9999)
            {   
                switch(p->F90)
                {
                    case 0:
                        break;
                    case 1:
                        advectWater_forBonnScheme(a,p,sweep);
                        break;
                    case 3:
                        advectPlane_NewWang(a,p,sweep);
                        break;
                }
                
            }
            else
            {
                switch(p->F90)
                {
                    case 0:
                        break;
                    case 1:
                        break;
                    case 3:
                        advectPlane_NewWang(a,p,sweep);
                        break;
                }
            }
        }
        
        pgc->start4(p,nx,1);
        pgc->start4(p,ny,1);
        pgc->start4(p,nz,1);
        pgc->start4(p,alpha,1);
        pgc->start4(p,V_w_m,1);
        pgc->start4(p,V_w_p,1);
        pgc->start4(p,phiS1,1);
        pgc->start4(p,vofS1,1);
        pgc->start4(p,phiS0,1);
        pgc->start4(p,vofS0,1);
        
        switch(p->F90)
        {
            case 0:
                break;
            case 1:
                transportVOF_Bonn(a,p,nSweep,sweep);
                break;
            case 2:
                transportVOF_NewWang(a,p,nSweep);
        }
        
        pgc->start4(p,phiS1,1);
        pgc->start4(p,vofS1,1);
        pgc->start4(p,phiS0,1);
        pgc->start4(p,vofS0,1);
        pgc->start4(p,vofS2,1);
        
        
        if(p->j_dir>0)
        {
             if(nSweep==2)
             {
                 pgc->start4(p,phiS2,1);
                 pgc->start4(p,vofS2,1);
                 LOOP
                 {
                     phistep(i,j,k)=phiS2(i,j,k);
                     vofstep(i,j,k)=vofS2(i,j,k);
                 }
                 pgc->start4(p,phistep,1);
                 pgc->start4(p,vofstep,1);
             }
        }
        else
        {
            if(nSweep==1)
            {
                pgc->start4(p,phiS1,1);
                pgc->start4(p,vofS1,1);
                LOOP
                {
                    phistep(i,j,k)=phiS1(i,j,k);
                    vofstep(i,j,k)=vofS1(i,j,k);
                }
                pgc->start4(p,phistep,1);
                pgc->start4(p,vofstep,1);
            }
        }
    }
    
    pgc->start4(p,phistep,1);
    pgc->start4(p,vofstep,1);
    
    LOOP
    {   
        if(p->F91<0.0)
        {
            if(vofstep(i,j,k)<0.0)
                vofstep(i,j,k)=0.0;
            if(vofstep(i,j,k)>1.0)
                vofstep(i,j,k)=1.0;
        }
        else
        {
            if(phistep(i,j,k)<-p->F91*p->psi || vofstep(i,j,k)<0.0)
                vofstep(i,j,k)=0.0;
            if(phistep(i,j,k)>p->F91*p->psi || vofstep(i,j,k)>1.0)
                vofstep(i,j,k)=1.0;
        }
            
        if(vofstep(i,j,k)>0.0001 && vofstep(i,j,k)<0.9999)
            reconstructPlane_alt(a,p,vofstep(i,j,k));
        else
        {
            nx(i,j,k)=2.0;
            ny(i,j,k)=0.0;
            nz(i,j,k)=2.0;
            alpha(i,j,k)=1E06;
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
            if((vofstep(i,j,k)>0.0001 && vofstep(i,j,k)<0.9999))
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
            else
                a->phi(i,j,k)=phistep(i,j,k);
        }
    
        pgc->start4(p,a->phi,1);
    
        reini_->start(a,p,a->phi,pgc,pflow);
    
        pgc->start4(p,a->phi,1);
    
    pupdate->start(p,a,pgc);
   
}

