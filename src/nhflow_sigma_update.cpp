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

#include"nhflow_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"

#define WLVL (fabs(WL(i,j))>0.00005?WL(i,j):1.0e20)

void nhflow_sigma::sigma_update(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    double wl,sigval;
    double bx,by,ex,ey;
    
    // calculate: Ex,Ey,Exx,Eyy
    // 3D
    if(p->i_dir==1 && p->j_dir==1)
    SLICELOOP4
    {
    d->Ex(i,j) = sx(d->eta);
    d->Ey(i,j) = sy(d->eta);
    
    d->Exx(i,j) = sxx(d->eta);
    d->Eyy(i,j) = syy(d->eta);
    }
    
    // 2D
    if(p->j_dir==0)
    SLICELOOP4
    {
    d->Ex(i,j) = sx(d->eta);    
    d->Exx(i,j) = sxx(d->eta);
    }
    
    SLICELOOP4
    if(p->wet[IJ]==0)
    {
    d->Ex(i,j)=0.0;
    d->Ey(i,j)=0.0;
    }
        
    pgc->gcsl_start4(p,d->Ex,1);
    pgc->gcsl_start4(p,d->Ey,1);
    
    // calculate: Bx,By,Bxx,Byy
    // 3D
    if(p->j_dir==1)
    SLICELOOP4
    {
    d->Bx(i,j) = sx(d->depth);
    d->By(i,j) = sy(d->depth);
    
    d->Bxx(i,j) = sxx(d->depth);
    d->Byy(i,j) = syy(d->depth);
    }

    // 2D
    if(p->j_dir==0 && p->A312!=1)
    SLICELOOP4
    {
    d->Bx(i,j) = sx(d->depth);    
    d->Bxx(i,j) = sxx(d->depth);
    }
    
    SLICELOOP4
    if(p->wet[IJ]==0)
    {
    d->Bx(i,j)=0.0;
    d->By(i,j)=0.0;
    }

    pgc->gcsl_start4(p,d->Bx,1);
    pgc->gcsl_start4(p,d->By,1);
    
    // -----------------------------------------------------
    
    // sigx
    FLOOP
    {
    if(p->wet[IJ]==0)
    p->sigx[FIJK] = 0.0;
    
    if(p->wet[IJ]==1)
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(d->Bx(i,j)/WLVL) - p->sig[FIJK]*(d->Ex(i,j)/WLVL);
    }
    
    // sigy
    FLOOP
    {
    if(p->wet[IJ]==0)
    p->sigy[FIJK] = 0.0;
    
    if(p->wet[IJ]==1)
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(d->By(i,j)/WLVL) - p->sig[FIJK]*(d->Ey(i,j)/WLVL);
    }    
    
    // sigz
    SLICELOOP4
    {
    if(p->wet[IJ]==0)
    p->sigz[IJ] = 0.0;
    
    if(p->wet[IJ]==1)
    p->sigz[IJ] = 1.0/WLVL;
    }

    // sigt
    FLOOP
    p->sigt[FIJK] = -(p->sig[FIJK]/WLVL)*d->detadt(i,j);

    // sigxx
    FLOOP
    if(p->wet[IJ]==1)
    {
        p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVL)*(d->Bxx(i,j) - pow(d->Bx(i,j),2.0)/WLVL) // xx
        
                      - (p->sig[FIJK]/WLVL)*(d->Exx(i,j) - pow(d->Ex(i,j),2.0)/WLVL)
                      
                      - (p->sigx[FIJK]/WLVL)*(d->Bx(i,j) + d->Ex(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(d->Bx(i,j)*d->Ex(i,j))
                      
                      
                      + ((1.0 - p->sig[FIJK])/WLVL)*(d->Byy(i,j) - pow(d->By(i,j),2.0)/WLVL) // yy
        
                      - (p->sig[FIJK]/WLVL)*(d->Eyy(i,j) - pow(d->Ey(i,j),2.0)/WLVL)
                      
                      - (p->sigy[FIJK]/WLVL)*(d->By(i,j) + d->Ey(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(d->By(i,j)*d->Ey(i,j));
    }
    
    FLOOP
    if(p->wet[IJ]==0)
    p->sigxx[FIJK]=0.0;
    
    // sig BC
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigx[FIJKm1] = p->sigx[FIJK];
            p->sigx[FIJKm2] = p->sigx[FIJK];
            p->sigx[FIJKm3] = p->sigx[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigx[FIJKp1] = p->sigx[FIJK];
            p->sigx[FIJKp2] = p->sigx[FIJK];
            p->sigx[FIJKp3] = p->sigx[FIJK];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigy[FIJKm1] = p->sigy[FIJK];
            p->sigy[FIJKm2] = p->sigy[FIJK];
            p->sigy[FIJKm3] = p->sigy[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigy[FIJKp1] = p->sigy[FIJK];
            p->sigy[FIJKp2] = p->sigy[FIJK];
            p->sigy[FIJKp3] = p->sigy[FIJK];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigxx[FIJKm1] = p->sigxx[FIJK];
            p->sigxx[FIJKm2] = p->sigxx[FIJK];
            p->sigxx[FIJKm3] = p->sigxx[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigxx[FIJKp1] = p->sigxx[FIJK];
            p->sigxx[FIJKp2] = p->sigxx[FIJK];
            p->sigxx[FIJKp3] = p->sigxx[FIJK];
        } 
    }
    
    
    FLOOP
    p->ZSN[FIJK] = p->ZN[KP]*d->WL(i,j) + d->bed(i,j);
    
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
    
    pgc->start7S(p,p->sigx,1);
    pgc->start7S(p,p->sigy,1);
    pgc->start7S(p,p->sigxx,1);
    pgc->start7S(p,p->sigt,1);
    pgc->start7S(p,p->ZSN,1);
    pgc->gcslparaxijk(p, p->sigz, 1);
}

