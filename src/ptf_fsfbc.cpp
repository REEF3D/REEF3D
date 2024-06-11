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

#include"ptf_fsfbc.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"onephase.h"
#include"fnpf_voiddisc.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"fnpf_cds6.h"
#include"fnpf_weno5.h"
#include"fnpf_weno7.h"
#include"sflow_bicgstab.h"
#include"hypre_struct2D.h"
#include"ptf_coastline.h"

ptf_fsfbc::ptf_fsfbc(lexer *p, fdm_ptf *a, ghostcell *pgc) : Fx(p),Fy(p),Fz(p),Ex(p),Ey(p),bx(p),by(p),Bx(p),By(p)
{    
    if(p->A311==0)
    pconvec = new fnpf_voiddisc(p);
    
    if(p->A311==2)
    pconvec = new fnpf_cds2(p);
    
    if(p->A311==3)
    pconvec = new fnpf_cds4(p);
    
    if(p->A311==4)
    pconvec = new fnpf_weno5(p);
    
    if(p->A311==6)
    pconvec = new fnpf_cds6(p);
    
    if(p->A311==7)
    pconvec = new fnpf_weno7(p);
    
    if(p->A350==1)
    psolv =  new sflow_bicgstab(p,pgc);
    
    pdx = new fnpf_cds2(p);

    a->wd_criterion=0.00005;
    
    if(p->A344==1)
    a->wd_criterion=p->A344_val;
    
    if(p->A345==1)
    a->wd_criterion=p->A345_val*p->DXM;
    
    pcoast = new ptf_coastline(p);
    
    dist3=0.0;
    
    if(p->A341>0.0 && p->j_dir==0)
    dist3=p->A341*(p->DXD);
    
    if(p->A341>0.0 && p->j_dir==1)
    dist3=p->A341*0.5*(p->DXD+p->DYD);
    
    if(p->A342>0.0)
    dist3=p->A342;
    
    dist4=1.0*dist3;
    
    expinverse = 1.0/(exp(1.0)-1.0);
    
    count_n=0;
    
}

ptf_fsfbc::~ptf_fsfbc()
{
}

void ptf_fsfbc::fsfdisc(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &eta, slice &Fifsf, field &Fi)
{   
    SLICELOOP4
        a->WL(i,j) = MAX(0.0, a->eta(i,j) + p->wd - a->bed(i,j));
    // 3D
    if(p->j_dir==1)
    {
        FILOOP4
        {
            ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);    
            jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    
            Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
            Fy(i,j) = pconvec->sy(p,Fifsf,jvel);
    
            Ex(i,j) = pconvec->sx(p,eta,ivel);
            Ey(i,j) = pconvec->sy(p,eta,jvel);
    
            Bx(i,j) = pdx->sx(p,a->depth,1.0);
            By(i,j) = pdx->sy(p,a->depth,1.0);
        }
    
    pgc->gcsl_start4(p,Ex,1);
    pgc->gcsl_start4(p,Ey,1);
    pgc->gcsl_start4(p,Bx,1);
    pgc->gcsl_start4(p,By,1);
    }
    
    // 2D
    if(p->j_dir==0)
    {
        FILOOP4
        {
            ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);      
    
            Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    
            Ex(i,j) = pconvec->sx(p,eta,ivel);
    
            Bx(i,j) = pdx->sx(p,a->depth,1.0);
        }
        pgc->gcsl_start4(p,Ex,1);
        pgc->gcsl_start4(p,Bx,1);
    }
}

void ptf_fsfbc::fsfdisc_ini(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    FILOOP4
    {
    Bx(i,j) = pdx->sx(p,a->depth,1.0);
    By(i,j) = pdx->sy(p,a->depth,1.0);
    }
    
    pgc->gcsl_start4(p,Bx,1);
    pgc->gcsl_start4(p,By,1);
    
    SLICELOOP4
        p->wet[IJ]=1;
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
}

void ptf_fsfbc::fsfwvel(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    FILOOP4
    {
    kvel = (a->Fi(i,j,k) - a->Fi(i,j,k-1))/(p->DZP[KP]);
    
    //if(p->A323<4)
    Fz(i,j) = pconvec->fz(p,a->Fi,kvel,kvel);
    
    if(p->wet[IJ]==0)
        Fz(i,j) = 0.0;
        
    coastline_fi(p,a,pgc,Fz);
    
    //if(p->A323==4)
    //Fz(i,j) = fz(p,a,a->Fi,Fifsf);
    
    //cout<<"Fz: "<<Fz(i,j)<<endl;
    //if(p->wet[IJ]==0)
    //c->Fz(i,j) = 0.0;
    }
}

void ptf_fsfbc::kfsfbc(lexer *p, fdm_ptf *a, ghostcell *pgc)
{
    SLICELOOP4
    {
    a->K(i,j) = - Fx(i,j)*Ex(i,j) - Fy(i,j)*Ey(i,j) 
    
                + Fz(i,j)*(1.0 + pow(Ex(i,j),2.0) + pow(Ey(i,j),2.0));
     
    if(p->wet[IJ]==0)
        a->K(i,j)=0.0;
    }
}

void ptf_fsfbc::dfsfbc(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &eta)
{
    SLICELOOP4
    {
    a->K(i,j) = - 0.5*pow(Fx(i,j),2.0) - 0.5*pow(Fy(i,j),2.0) 
    
                + 0.5*pow(Fz(i,j),2.0)*(1.0 + pow(Ex(i,j),2.0) + pow(Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j);
    
    if(p->wet[IJ]==0)
        a->K(i,j)=0.0;
    }
}


