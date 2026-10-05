/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"nhflow_breaking.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_breaking::nhflow_breaking(lexer* p, fdm_nhf *d, ghostcell *pgc) : bx(p), by(p), bd(p), brkflag(p), count_n(-1)
{
    SLICELOOP4
    d->breaking(i,j)=0;
}

nhflow_breaking::~nhflow_breaking()
{
}

void nhflow_breaking::breaking(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    if(p->A550==1)
    breaking_baquet(p, d, pgc, eta, eta_n, WL, alpha);
}

void nhflow_breaking::breaking_baquet(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    // mesh refinement (nhflow_amr, brk_defer): the three parts run on all grids in turn, the
    // patches exchange the detection and the flags of their sibling cells in between
    if(brk_defer)
    {
    brk_eta = &eta;
    brk_eta_n = &eta_n;
    brk_WL = &WL;
    brk_alpha = alpha;
    return;
    }

    breaking_detect(p,d,pgc,eta,eta_n,WL,alpha);
    breaking_flags(p,d,pgc);
    breaking_vb(p,d,pgc,WL);
}

void nhflow_breaking::breaking_part(lexer *p, fdm_nhf *d, ghostcell *pgc, int part)
{
    if(part==0)
    breaking_detect(p,d,pgc,*brk_eta,*brk_eta_n,*brk_WL,brk_alpha);
    if(part==1)
    breaking_flags(p,d,pgc);
    if(part==2)
    breaking_vb(p,d,pgc,*brk_WL);
}

// detection: steepness (bx, by: 10 / 20) and depth induced (bd: 2)
void nhflow_breaking::breaking_detect(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    if(p->count>count_n)
    {
    SLICELOOP4
    {
    brkflag(i,j)=0;
    d->breaking(i,j)=0;
    }
    
    count_n=p->count;
    }
    
    SLICELOOP4
    {
    bx(i,j)=0;
    by(i,j)=0;
    bd(i,j)=0;
    }
    
    pgc->gcsl_start4int(p,brkflag,50);
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    pgc->gcsl_start4int(p,bd,50);
    
    // steepness induced breaking
    if((p->A551==2 || p->A551==3) && p->count>1)
    {
    SLICELOOP4
    {
            // x
            if(d->Ex(i,j)   < -p->A555)
            {
                bx(i,j) = 10;
            }
            
            if(d->Ex(i,j)   > p->A555)
            {
                bx(i,j) = 20;

            }
            
            // y
            if(p->j_dir==1)
            if( d->Ey(i,j)   < -p->A555)
            {
                by(i,j) = 10;
            }
            
            if(p->j_dir==1)
            if( d->Ey(i,j)   > p->A555)
            {
                by(i,j) = 20;
            }
            
    }
    
    // mesh refinement: no detection in the filled cells around a patch (their surface comes
    // from the coarser grid and their outer cells see one-sided fluxes); their breaking comes
    // from the source grid below
    if(p->amrvb!=nullptr)
    SLICELOOP4
    if(p->amrvb[IJ]>=0.0)
    {
    bx(i,j)=0;
    by(i,j)=0;
    }

    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    }
    
    // depth induced breaking
    if((p->A551==1 || p->A551==3) && p->count>1)
    {
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if((eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A554*sqrt(9.81*d->WL(i,j)))
        bd(i,j)=2;
    }
    
    if(p->amrvb!=nullptr)
    SLICELOOP4
    if(p->amrvb[IJ]>=0.0)
    bd(i,j)=0;

    pgc->gcsl_start4int(p,bd,50);
    }
}

// the breaking flags from the detection
void nhflow_breaking::breaking_flags(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if((p->A551==2 || p->A551==3) && p->count>1)
    {
    SLICELOOP4
    if(bx(i,j)>0 || by(i,j)>0)
    {
    brkflag(i,j)=1;
    }
    
    SLICELOOP4
    {
    // x
    if(bx(i,j)==10 || bx(i-1,j)==10 || bx(i+1,j)==10 || bx(i+2,j)==10)
    brkflag(i,j)=1;
    
    if(bx(i,j)==20 || bx(i-1,j)==20 || bx(i-2,j)==20 || bx(i+1,j)==20)
    brkflag(i,j)=1;
    
    // y
    if(by(i,j)==10 || by(i,j-1)==10 || by(i,j+1)==10 || by(i,j+2)==10)
    brkflag(i,j)=1;
    
    if(by(i,j)==20 || by(i,j-1)==20 || by(i,j-2)==20 || by(i,j+1)==20)
    brkflag(i,j)=1;
        
    }
    }
    
    if((p->A551==1 || p->A551==3) && p->count>1)
    SLICELOOP4
    {
    // x
    if(bd(i,j)==2 || bd(i-1,j)==2 || bd(i-2,j)==2 || bd(i+1,j)==2 || bd(i+2,j)==2)
    brkflag(i,j)=2;
    
    // y
    if(p->j_dir==1)
    if(bd(i,j)==2 || bd(i,j-1)==2 || bd(i,j-2)==2 || bd(i,j+1)==2 || bd(i,j+2)==2)
    brkflag(i,j)=2;    
    }
    
    // mesh refinement: a filled cell breaks where its source cell carries the full breaking
    // viscosity
    if(p->amrvb!=nullptr)
    SLICELOOP4
    if(p->amrvb[IJ]>=p->A557 && brkflag(i,j)==0)
    brkflag(i,j)=1;
}

// the breaking viscosity from the flags
void nhflow_breaking::breaking_vb(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    // ------------------------
    // fill flag
    // ------------------------
    /*if(p->j_dir==0)
    SLICELOOP4
    {   
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0)
        d->breaking(i,j) = 1;
            
        if(brkflag(i,j)==0 && 
        (brkflag(i-1,j)>0 || brkflag(i+1,j)>0 || brkflag(i-2,j)>0 || brkflag(i+2,j)>0))
        d->breaking(i,j) = 1;
    }

    if(p->j_dir==1)
    SLICELOOP4
    {   
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0 || brkflag(i,j-1)>0 || brkflag(i,j+1)>0)
        d->breaking(i,j) = 1;
            
        if(brkflag(i,j)==0 && 
        (brkflag(i-1,j-1)>0 || brkflag(i-1,j+1)>0 || brkflag(i+1,j-1)>0 || brkflag(i+1,j+1)>0
        || brkflag(i-2,j)>0 || brkflag(i+2,j)>0 || brkflag(i,j-2)>0 || brkflag(i,j+2)>0))
        d->breaking(i,j) = 1;
    }
    */
    
    
    
    
    // ------------------------
    // fill breaking viscosity
    // ------------------------
    
    SLICELOOP4
    d->vb(i,j) = 0.0;
        
    if(p->j_dir==0)
    SLICELOOP4
    {   
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0)
        d->vb(i,j) = p->A557;
            
        if(brkflag(i,j)==0 && 
        (brkflag(i-1,j)>0 || brkflag(i+1,j)>0 || brkflag(i-2,j)>0 || brkflag(i+2,j)>0))
        d->vb(i,j) = 0.5*p->A557;
    }

    if(p->j_dir==1)
    SLICELOOP4
    {   
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0 || brkflag(i,j-1)>0 || brkflag(i,j+1)>0)
        d->vb(i,j) = p->A557;
            
        if(brkflag(i,j)==0 && 
        ( brkflag(i-1,j-1)>0 || brkflag(i-1,j+1)>0 || brkflag(i+1,j-1)>0 || brkflag(i+1,j+1)>0
        || brkflag(i-2,j)>0 || brkflag(i+2,j)>0 || brkflag(i,j-2)>0 || brkflag(i,j+2)>0))
        d->vb(i,j) = 0.5*p->A557;
    }
    
        
        // additional breaking filter
        // shallow
        if(p->A552==1)
        SLICELOOP4
        if(brkflag(i,j)==2)
        {
         filter(p,d,pgc,WL);
        }   
        
        // deep
        if(p->A552==2)
        SLICELOOP4
        if(brkflag(i,j)==1)
        {
         filter(p,d,pgc,WL);
        }   
        
        // all
        if(p->A552==3)
        SLICELOOP4
        if(brkflag(i,j)>=1)
        {
         filter(p,d,pgc,WL);
        }   
        
        // coastline filter
        /*SLICELOOP4
        {
            
            if(d->coastline(i,j)>=0.0)
            {
                db = d->coastline(i,j);
                
                if(db<dist3)
                {
                filter(p,d,pgc,eta);
                filter(p,d,pgc,Fifsf);
                }
            }
        }*/
        
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4int(p,d->breaking,50);
    
    SLICELOOP4
    d->breaklog(i,j)=0;
    
    // breaklog
    int count=0;
    
    SLICELOOP4
    if(d->breaking(i,j)>0)
    {
    d->breaklog(i,j)=1;
    ++count;
    }

    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}
