/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_sg_fsfbc_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"onephase.h"
#include"fnpf_voiddisc.h"
#include"fnpf_cds2_wd.h"
#include"fnpf_cds4_wd.h"
#include"fnpf_cds6_wd.h"
#include"fnpf_weno.h"
#include"fnpf_weno_wd.h"
#include"fnpf_weno7.h"
#include"fnpf_wenoflux.h"
#include"fnpf_ddx_cds2_wd.h"
#include"fnpf_ddx_cds4_wd.h"
#include"fnpf_sg_coastline.h"
#include"sflow_bicgstab.h"

fnpf_sg_fsfbc_wd::fnpf_sg_fsfbc_wd(lexer *p, fdm_fnpf *c, ghostcell *pgc) : bx(p),by(p)
{    
    if(p->A311==0)
    pconvec = pconeta = new fnpf_voiddisc(p);
    
    if(p->A311==1)
    pconvec = pconeta = new fnpf_cds2_wd(p,c);
    
    if(p->A311==2)
    pconvec = pconeta = new fnpf_cds4_wd(p);
    
    if(p->A311==4  || p->A311==5)
    {
    pconvec = new fnpf_weno_wd(p,c);
    pconeta = new fnpf_weno_wd(p,c);
    }
    //pdf = new fnpf_wenoflux(p);

    if(p->A311==6)
    pconvec = pconeta = new fnpf_cds6_wd(p);
    
    if(p->A311==7)
    {
    pconvec = new fnpf_weno7(p);
    pconeta = new fnpf_weno7(p);
    }
    

    // ---
    if(p->A312==2)
    {
    pddx = new fnpf_ddx_cds2_wd(p,c);
    pdx = new fnpf_cds2_wd(p,c);
    }
    
    if(p->A312==3)
    {
    pddx = new fnpf_ddx_cds4_wd(p);
    pdx = new fnpf_cds4_wd(p);
    }
    
    
    FFILOOP4
    {
    c->Fy(i,j) = 0.0;
    c->Ey(i,j) = 0.0;
    c->Hy(i,j) = 0.0;
    c->Eyy(i,j) = 0.0;
    }
    
    
    c->wd_criterion=0.00005;
    
    if(p->A344==1)
    c->wd_criterion=p->A244_val;
    
    if(p->A345==1)
    c->wd_criterion=p->A245_val*p->dx;
    
    pcoast = new fnpf_sg_coastline(p);
    
    dist3=0.0;
    
    if(p->A341>0.0)
    dist3=p->A341*p->DXM;
    
    if(p->A342>0.0)
    dist3=p->A342;
    
    expinverse = 1.0/(exp(1.0)-1.0);
    
    if(p->A350==1)
    psolv = new sflow_bicgstab(p,pgc);
    
    count_n=0;
}

fnpf_sg_fsfbc_wd::~fnpf_sg_fsfbc_wd()
{
}

void fnpf_sg_fsfbc_wd::fsfdisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    c->WL(i,j) = MAX(0.0,c->eta(i,j) + p->wd - c->bed(i,j));
    
    pgc->gcsl_start4(p,c->WL,50);
    
    // 3D
    if(p->i_dir==1 && p->j_dir==1)
    FFILOOP4
    {
    ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);    
    jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    
    c->Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    c->Fy(i,j) = pconvec->sy(p,Fifsf,jvel);
    
    c->Ex(i,j) = pconeta->sx(p,eta,ivel);
    c->Ey(i,j) = pconeta->sy(p,eta,jvel);
    
    c->Exx(i,j) = pddx->sxx(p,eta);
    c->Eyy(i,j) = pddx->syy(p,eta);
    }
    
    // 2D
    if(p->i_dir==1 && p->j_dir==0)
    FFILOOP4
    {
    ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);    
    
    c->Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    c->Ex(i,j) = pconeta->sx(p,eta,ivel);
    
    c->Exx(i,j) = pddx->sxx(p,eta);
    }
}

void fnpf_sg_fsfbc_wd::fsfdisc_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    FFILOOP4
    {
    c->Bx(i,j) = pdx->sx(p,c->depth,1.0);
    c->By(i,j) = pdx->sy(p,c->depth,1.0);
    
    c->Bxx(i,j) = pddx->sxx(p,c->depth);
    c->Byy(i,j) = pddx->syy(p,c->depth);
    }
    
    pgc->gcsl_start4(p,c->Bx,1);
    pgc->gcsl_start4(p,c->By,1);
    
    SLICELOOP4
    c->wet(i,j)=1;
    
    pgc->gcsl_start4int(p,c->wet,50);
}

void fnpf_sg_fsfbc_wd::fsfwvel(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    FFILOOP4
    {
    c->Fz(i,j) = p->sigz[IJ]*pconvec->sz(p,c->Fi);

    if(c->wet(i,j)==0)
    c->Fz(i,j) = 0.0;
    }
}

void fnpf_sg_fsfbc_wd::kfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    double dEdF_x,dEdF_y;
    
    SLICELOOP4
    {
    dEdF_x = c->Fx(i,j)*c->Ex(i,j);
    dEdF_y = c->Fy(i,j)*c->Ey(i,j);

    c->K(i,j) =  - dEdF_x - dEdF_y
    
                 + c->Fz(i,j)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0));
                 
          
    if(c->wet(i-1,j)==0 || c->wet(i+1,j)==0 || c->wet(i,j-1)==0 || c->wet(i,j+1)==0 
  || c->wet(i-1,j-1)==0 || c->wet(i+1,j-1)==0 || c->wet(i-1,j+1)==0 || c->wet(i+1,j+1)==0)
      c->K(i,j) = c->Fz(i,j);
    }
}

void fnpf_sg_fsfbc_wd::dfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta)
{
    double dFdF_x,dFdF_y;
    
    SLICELOOP4
    {
    dFdF_x = c->Fx(i,j)*c->Fx(i,j);
    dFdF_y = c->Fy(i,j)*c->Fy(i,j);
    
    c->K(i,j) =  - 0.5*dFdF_x - 0.5*dFdF_y
    
                 + 0.5*pow(c->Fz(i,j),2.0)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j);
                 
            
    if(c->wet(i-1,j)==0 || c->wet(i+1,j)==0 || c->wet(i,j-1)==0 || c->wet(i,j+1)==0 
  || c->wet(i-1,j-1)==0 || c->wet(i+1,j-1)==0 || c->wet(i-1,j+1)==0 || c->wet(i+1,j+1)==0)
      c->K(i,j) =  - fabs(p->W22)*eta(i,j);
    }
}

