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

#include"fnpf_fsfbc.h"
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
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"fnpf_cds6.h"
#include"fnpf_weno3.h"
#include"fnpf_weno5.h"
#include"fnpf_weno7.h"
#include"fnpf_weno5_wd.h"
#include"fnpf_wenoflux.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"
#include"sflow_bicgstab.h"
#include"hypre_struct2D.h"

fnpf_fsfbc::fnpf_fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc) : bx(p),by(p),eps(1.0e-6)
{    
    if(p->A311==0)
    {
    pconvec = new fnpf_voiddisc(p);
    pconeta = new fnpf_voiddisc(p);
    }
    
    if(p->A311==1)
    {
    pconvec = new fnpf_cds2(p);
    pconeta = new fnpf_cds2(p);
    }
    
    if(p->A311==2)
    {
    pconvec = new fnpf_cds4(p);
    pconeta = new fnpf_cds4(p);
    }
    
    if(p->A311==3)
    {
    pconvec = new fnpf_weno3(p);
    pconeta = new fnpf_weno3(p);
    }
    
    if(p->A311==4 || p->A311==5)
    {
    pconvec = new fnpf_weno5(p);
    pconeta = new fnpf_weno5(p);
    }
    
    pdf = new fnpf_wenoflux(p);

    if(p->A311==6)
    {
    pconvec = new fnpf_cds6(p);
    pconeta = new fnpf_cds6(p);
    }
    
    if(p->A311==7)
    {
    pconvec = new fnpf_weno7(p);
    pconeta = new fnpf_weno7(p);
    }
    

    // ---
    if(p->A312==2)
    {
    pddx = new fnpf_ddx_cds2(p);
    pdx = new fnpf_cds2(p);
    }
    
    if(p->A312==3)
    {
    pddx = new fnpf_ddx_cds4(p);
    pdx = new fnpf_cds4(p);
    }
    
    
    FFILOOP4
    {
    c->Fy(i,j) = 0.0;
    c->Ey(i,j) = 0.0;
    c->Hy(i,j) = 0.0;
    c->Eyy(i,j) = 0.0;
    c->By(i,j) = 0.0;
    c->Byy(i,j) = 0.0;
    }
    
    
    c->wd_criterion=0.00005;
    
    if(p->A344==1)
    c->wd_criterion=p->A244_val;
    
    if(p->A345==1)
    c->wd_criterion=p->A245_val*p->DXM;
    
    if(p->A350==1)
    psolv =  new sflow_bicgstab(p,pgc);
    
    
    p->Iarray(temp,p->imax*p->jmax);
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
}

fnpf_fsfbc::~fnpf_fsfbc()
{
}

void fnpf_fsfbc::fsfdisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    c->WL(i,j) = MAX(0.0, c->eta(i,j) + p->wd - c->bed(i,j));

    
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
    
    pgc->gcsl_start4(p,c->Ex,1);
    pgc->gcsl_start4(p,c->Ey,1);
}

void fnpf_fsfbc::fsfdisc_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // 3D
    if(p->i_dir==1 && p->j_dir==1)
    FFILOOP4
    {
    c->Bx(i,j) = pdx->sx(p,c->depth,1.0);
    c->By(i,j) = pdx->sy(p,c->depth,1.0);
    
    c->Bxx(i,j) = pddx->sxx(p,c->depth);
    c->Byy(i,j) = pddx->syy(p,c->depth);
    }
    
    // 2D
    if(p->i_dir==1 && p->j_dir==0)
    FFILOOP4
    {
    c->Bx(i,j) = pdx->sx(p,c->depth,1.0);    
    c->Bxx(i,j) = pddx->sxx(p,c->depth);
    }
    
    pgc->gcsl_start4(p,c->Bx,1);
    pgc->gcsl_start4(p,c->By,1);
    
    SLICELOOP4
    p->wet[IJ]=1;
    
    pgc->gcsl_start4Vint(p,p->wet,50);
}

void fnpf_fsfbc::fsfwvel(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    FFILOOP4
    {
    c->Fz(i,j) = p->sigz[IJ]*pconvec->sz(p,c->Fi);
    
    if(p->wet[IJ]==0)
    c->Fz(i,j) = 0.0;
    }
}

void fnpf_fsfbc::kfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    SLICELOOP4
    c->K(i,j) =  - c->Fx(i,j)*c->Ex(i,j) - c->Fy(i,j)*c->Ey(i,j)
    
                 + c->Fz(i,j)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0));
}

void fnpf_fsfbc::dfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta)
{ 
    SLICELOOP4
    c->K(i,j) =  - 0.5*c->Fx(i,j)*c->Fx(i,j) - 0.5*c->Fy(i,j)*c->Fy(i,j)
    
                 + 0.5*pow(c->Fz(i,j),2.0)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j);
}


void fnpf_fsfbc::wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf) 
{   
    
    SLICELOOP4
    p->wet[IJ]=1;
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    /*

    SLICELOOP4
    c->WL(i,j) = eta(i,j) + p->wd - c->bed(i,j);

    
    pgc->gcsl_start4(p,c->WL,50);
      
    
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
     
    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
            if(p->wet[Ip1J]==1 && eta(i,j)<eta(i+1,j) && c->WL(i+1,j)>c->wd_criterion+eps)
            temp[IJ]=1;
            
            if(p->wet[Im1J]==1 && eta(i,j)<eta(i-1,j) && c->WL(i-1,j)>c->wd_criterion+eps)
            temp[IJ]=1;
            
            if(p->wet[IJp1]==1 && eta(i,j)<eta(i,j+1) && c->WL(i,j+1)>c->wd_criterion+eps && p->j_dir==1)
            temp[IJ]=1;
            
            if(p->wet[IJm1]==1 && eta(i,j)<eta(i,j-1) && c->WL(i,j-1)>c->wd_criterion+eps && p->j_dir==1)
            temp[IJ]=1;
        }
        
        else              
        if(c->WL(i,j)<=c->wd_criterion)
        {
        temp[IJ]=0;
        eta(i,j) = c->wd_criterion - c->depth(i,j);
        c->WL(i,j) = eta(i,j) + c->depth(i,j);
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];
    

    pgc->gcsl_start4Vint(p,p->wet,50);
    pgc->gcsl_start4(p,eta,gcval_eta);
    pgc->gcsl_start4(p,c->WL,gcval_eta);*/
      
    SLICELOOP4
    c->test2D(i,j) = double (p->wet[IJ]);
}

void fnpf_fsfbc::depthcheck(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf) 
{ 
}

void fnpf_fsfbc::coastline_eta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f) 
{  
}

void fnpf_fsfbc::coastline_fi(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f) 
{   
}
