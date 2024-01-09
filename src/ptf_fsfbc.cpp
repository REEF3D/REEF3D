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

#include"ptf_fsfbc.h"
#include"lexer.h"
#include"fdm.h"
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

ptf_fsfbc::ptf_fsfbc(lexer *p, fdm *a, ghostcell *pgc) : Fx(p),Fy(p),Fz(p),Ex(p),Ey(p)
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
    
}

ptf_fsfbc::~ptf_fsfbc()
{
}

void ptf_fsfbc::fsfdisc(lexer *p, fdm *a, ghostcell *pgc, slice &eta, slice &Fifsf, field &Fi)
{
    // 3D
    if(p->j_dir==1)
    FILOOP4
    {
    ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);    
    jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    
    Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    Fy(i,j) = pconvec->sy(p,Fifsf,jvel);
    
    Ex(i,j) = pconvec->sx(p,eta,ivel);
    Ey(i,j) = pconvec->sy(p,eta,jvel);
    }
    
    // 2D
    if(p->j_dir==0)
    FILOOP4
    {
    ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);      
    
    Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    
    Ex(i,j) = pconvec->sx(p,eta,ivel);
    }
    
    pgc->gcsl_start4(p,Ex,1);
    pgc->gcsl_start4(p,Ey,1);
}

void ptf_fsfbc::fsfwvel(lexer *p, fdm *a, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    FILOOP4
    {
    kvel = (a->Fi(i,j,k) - a->Fi(i,j,k-1))/(p->DZP[KP]);
    
    //if(p->A323<4)
    Fz(i,j) = pconvec->fz(p,a->Fi,kvel,kvel);
    
    //if(p->A323==4)
    //Fz(i,j) = fz(p,a,a->Fi,Fifsf);
    
    //cout<<"Fz: "<<Fz(i,j)<<endl;
    //if(p->wet[IJ]==0)
    //c->Fz(i,j) = 0.0;
    }
}

void ptf_fsfbc::kfsfbc(lexer *p, fdm *a, ghostcell *pgc)
{
    SLICELOOP4
    a->K(i,j) = - Fx(i,j)*Ex(i,j) - Fy(i,j)*Ey(i,j) 
    
                + Fz(i,j)*(1.0 + pow(Ex(i,j),2.0) + pow(Ey(i,j),2.0));
}

void ptf_fsfbc::dfsfbc(lexer *p, fdm *a, ghostcell *pgc, slice &eta)
{
    SLICELOOP4
    a->K(i,j) = - 0.5*pow(Fx(i,j),2.0) - 0.5*pow(Fy(i,j),2.0) 
    
                + 0.5*pow(Fz(i,j),2.0)*(1.0 + pow(Ex(i,j),2.0) + pow(Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j);
}


