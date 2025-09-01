/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include"picard_f.h"


VOF_PLIC::VOF_PLIC
(
    lexer* p,
    fdm *a,
    ghostcell* pgc,
    heat *pheat
):gradient(p),norm_vec(p),alpha(p),nx(p),ny(p),nz(p),vof1(p),vof2(p),vof3(p),phival(p),V_w_p(p),V_w_m(p),alphastore(p),phistep(p),phiS0(p),phiS1(p),phiS2(p),vofstep(p),vofS0(p),vofS1(p),vofS2(p),
    phiaux(p),Vx_p(p),Vx_m(p),Vz_p(p),Vz_m(p),Vn_p(p),Vn_m(p),F_x(p),F_z(p),F_n(p),F_new(p),Flux_x(p),Flux_z(p),Crossflux_xz(p),Crossflux_zx(p),curv(p),compressvol(p),VoF(p),vof_rk1(p),vof_rk2(p)
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

    preini = new reini_RK3(p,1);
    
    ppicard = new picard_f(p);
    
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
    
    a_thres=p->F93;
    w_thres=p->F94;
    corr_thres=p->F95;
    gcval_vof=1;
    gcval_ro=1;
    gcval_visc=1;
    gcval_phi=gcval_frac;
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
    pupdate->start(p,a,pgc,a->u,a->v,a->w);
}

void VOF_PLIC::start(fdm* a,lexer* p, convection* pconvec,solver* psolv, ghostcell* pgc,ioflow* pflow, reini* preini, particle_corr* ppart, field &ls)
{	
    
//********************************************************
//Step 1
//********************************************************
    // get vectorized face density from density_f
    
    //-------------------------------------------
    // FSF
    LOOP
    {
	a->L(i,j,k)=0.0;
    VoF(i,j,k)=a->vof(i,j,k);
    }
/*
	RKcalcL(a,p,pgc,a->u,a->v,a->w);
	
	LOOP
    {
        vof_rk1(i,j,k) = VoF(i,j,k) + a->L(i,j,k);
        
        if(vof_rk1(i,j,k)<0.0)
            vof_rk1(i,j,k)=0.0;
        if(vof_rk1(i,j,k)>1.0)
            vof_rk1(i,j,k)=1.0;
    }
    
    pgc->start4(p,vof_rk1,gcval_vof);
    updatePlaneData(p,a,pgc,vof_rk1);
	pflow->vof_relax(p,a,pgc,vof_rk1);
	pgc->start4(p,vof_rk1,gcval_vof);
    
    LOOP
    {
     a->vof(i,j,k) = vof_rk1(i,j,k);
     a->L(i,j,k)=0.0;
    }
    pgc->start4(p,a->vof,gcval_vof); 
    
    //!no update yet -> update after diffusion!
    
    if(p->F92==1)
    {
        RK_redistance(a,p,pgc);
        pgc->start4(p,a->phi,gcval_phi);
        
        p->F44=3;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
    }
    
    //-------------------------------------------
    
    updatePhasemarkersCompression(p,a,pgc,vof_rk1);
    pgc->start4(p,vof_rk1,gcval_vof);
//********************************************************
//Step 2
//********************************************************

    //-------------------------------------------
    // FSF
    
   RKcalcL(a,p,pgc,a->u,a->v,a->w);
	
	LOOP
    {
        vof_rk2(i,j,k) = 0.75*VoF(i,j,k) + 0.25*vof_rk1(i,j,k)+0.25*a->L(i,j,k);
        
        if(vof_rk2(i,j,k)<0.0)
            vof_rk2(i,j,k)=0.0;
        if(vof_rk2(i,j,k)>1.0)
            vof_rk2(i,j,k)=1.0;
    }
    
    updatePlaneData(p,a,pgc,vof_rk2);
	pflow->vof_relax(p,a,pgc,vof_rk2);
    pgc->start4(p,vof_rk2,gcval_vof);
    
    LOOP
    {
        a->vof(i,j,k) = vof_rk2(i,j,k);
        a->L(i,j,k)=0.0;
    }
    pgc->start4(p,a->vof,gcval_vof);
    
     if(p->F92==1)
    {
        RK_redistance(a,p,pgc);
        pgc->start4(p,a->phi,gcval_phi);
        
        p->F44=3;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
    }
    
    updatePhasemarkersCompression(p,a,pgc,vof_rk2);
    pgc->start4(p,vof_rk2,gcval_vof);
    
//********************************************************
//Step 3
//********************************************************
    //-------------------------------------------
    // FSF
    */
    RKcalcL(a,p,pgc,a->u,a->v,a->w);
	
	LOOP
    {
        a->vof(i,j,k) = VoF(i,j,k) + a->L(i,j,k);
        
        if(a->vof(i,j,k)<0.0)
            a->vof(i,j,k)=0.0;
        if(a->vof(i,j,k)>1.0)
            a->vof(i,j,k)=1.0;
    }
    
    //updatePlaneData(p,a,pgc,a->vof);
    pflow->vof_relax(p,a,pgc,a->vof);
    pgc->start4(p,a->vof,gcval_vof);
    
    LOOP
        a->L(i,j,k)=0.0;
    
    if(p->F92==1)
    {
        RK_redistance(a,p,pgc);
        pgc->start4(p,a->phi,gcval_phi);
        
        p->F44=4;
        preini->start(a,p,a->phi, pgc, pflow);
        ppicard->correct_ls(p,a,pgc,a->phi);
    }
    
    //-------------------------------------------
        
    updatePhasemarkersCorrection(p,a,pgc,a->vof);
    pgc->start4(p,a->vof,gcval_vof);
    //updatePlaneData(p,a,pgc,a->vof);
    
    if(p->F92==3)
        calculateSubFractions(p,a,pgc,a->vof);
    pupdate->start(p,a,pgc,a->u,a->v,a->w);
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    
    LOOP
    {
        if(a->vof(i,j,k)>p->F94)
            a->phi(i,j,k)=1.0;
        else if(a->vof(i,j,k)<p->F93)
            a->phi(i,j,k)=-1.0;
        else
            a->phi(i,j,k)=(a->vof(i,j,k)-0.5)*p->DZN[KP];
    }
    pgc->start4(p,a->phi,1);
    
    
    double vofchecksum;
    vofchecksum=0.0;
    LOOP
        vofchecksum+=a->vof(i,j,k)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
    vofchecksum=pgc->globalsum(vofchecksum);
    cout<<"Total water volume:"<<vofchecksum<<endl;
}
