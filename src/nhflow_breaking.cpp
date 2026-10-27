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
#include"ioflow.h"
#include<vector>
#include"breaking_crest3d.h"

nhflow_breaking::nhflow_breaking(lexer* p, fdm_nhf *d, ghostcell *pgc) : bx(p), by(p), bd(p), brkflag(p), brkzone(p),
                                                                            brkage(p), brkage_n(p), brkvb(p), bt_us(p), bt_vs(p), bt_vbs(p), brkseed(p),
                                                                            brkkseed(p), brkk(p), brkk_n(p), brkr(p), brkr_n(p), count_n(-1)
{
    SLICELOOP4
    {
    d->breaking(i,j)=0;
    brkflag(i,j)=0;
    brkzone(i,j)=0;
    brkage(i,j)=0.0;
    brkage_n(i,j)=0.0;
    brkseed(i,j)=0;
    brkkseed(i,j)=0.0;
    brkk(i,j)=0.0;
    brkk_n(i,j)=0.0;
    brkr(i,j)=0.0;
    brkr_n(i,j)=0.0;
    }
}

nhflow_breaking::~nhflow_breaking()
{
    delete bt;
}

void nhflow_breaking::breaking(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, slice &WL, double alpha)
{
    if(p->A550==1)
    breaking_baquet(p, d, pgc, eta, eta_n, WL, alpha);
    
    if(p->A550==3)
    breaking_kennedy(p, d, pgc, WL);
    
    if(p->A550==4)
    {
    breaking_barthelemy(p, d, pgc, eta_n);
    breaking_kennedy(p, d, pgc, WL);
    }
}

// A 550 3: eddy viscosity of Kennedy, Chen, Kirby & Dalrymple (2000), J. Waterway Port Coastal Ocean Eng. 126(1), 39-47,
// as in FUNWAVE, with deta/dt of the stage (continuity) and the total depth h of the stage water level:
//
//   nu_b = B delta_b^2 h deta/dt,   delta_b = A 504 (1.2)
//   B = 1 for deta/dt >= 2 eta_t*,  deta/dt / eta_t* - 1 for eta_t* < deta/dt < 2 eta_t*,  0 below
//   eta_t* = eta_t^I + min(t_b/T*, 1) (eta_t^F - eta_t^I),  eta_t^I = A 554 sqrt(g h), eta_t^F = A 501 sqrt(g h),
//   T* = A 503 sqrt(h/g)
//
// t_b is the age of the breaking event in the column: a column with deta/dt > eta_t* takes the largest age of
// itself and its neighbours at the end of the last time step plus dt, so the event moves with the front and its
// threshold relaxes from eta_t^I to eta_t^F while it travels; otherwise the age is 0. The viscosity varies
// continuously with deta/dt (B), so there is no switch of the flow equations and no on/off of a constant viscosity.
//
// nu_b is spread (maximum of the neighbours) over a band of A 502 h around the breaking columns, so that it covers
// the crest behind the front: on the front face alone it smooths the front, which takes away the dissipation of the
// shock-capturing scheme, while the non-hydrostatic pressure keeps the crest, and the crest peaks (Ting & Kirby).
// It acts in all layers on all faces of the implicit diffusion (A 512 2, set by A 550 3); with the horizontal faces
// only, the inner surf zone of Ting & Kirby depended on the number of layers (4: H/h 0.48, 10: 0.59 and noisy).
// No breaking in the relaxation zones (breaking_zones).
void nhflow_breaking::breaking_kennedy(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    const double g = fabs(p->W22);
    
    // new time step: the ages (and with A 550 4 the wavenumbers) at the end of the last step
    if(p->count>count_n)
    {
    SLICELOOP4
    {
    brkage_n(i,j) = brkage(i,j);
    brkk_n(i,j) = brkk(i,j);
    brkr_n(i,j) = brkr(i,j);
    }
    
    pgc->gcsl_start4(p,brkage_n,1);
    pgc->gcsl_start4(p,brkk_n,1);
    pgc->gcsl_start4(p,brkr_n,1);
    
    count_n=p->count;
    }
    
    const bool seeded = (p->A550==4);
    int count=0;
    
    SLICELOOP4
    {
    brkage(i,j) = 0.0;
    brkk(i,j) = 0.0;
    brkr(i,j) = 0.0;
    brkflag(i,j) = 0;
    d->vb(i,j) = 0.0;
    
        if(p->wet[IJ]==1 && brkzone(i,j)==0 && p->count>1)
        {
        const double h = MAX(WL(i,j),0.0);
        const double et = d->detadt(i,j);
        
        // the largest age of the column and its neighbours, and the wavenumber and initial threshold of that event (A 550 4)
        double age = brkage_n(i,j);
        double kw = brkk_n(i,j);
        double r0 = brkr_n(i,j);
        
        auto take = [&](int ii, int jj) { if(brkage_n(ii,jj)>age) { age = brkage_n(ii,jj); kw = brkk_n(ii,jj); r0 = brkr_n(ii,jj); } };
        take(i-1,j);
        take(i+1,j);
        
        if(p->j_dir==1)
        {
        take(i,j-1);
        take(i,j+1);
        }
        
        // A 550 4: a seeded column starts or continues an event with the wavenumber of the breaking crest; without
        // an event and without seed no breaking
        bool seedc = false;
        
        if(seeded)
        {
            if(brkseed(i,j)==1)
            {
                if(age<=0.0 || kw<=0.0)
                kw = brkkseed(i,j);
                
                seedc = true;
            }
            else
            if(age<=0.0)
            continue;
        }
        
        // length of the wave: h (A 550 3), tanh(k h)/k (A 550 4, deep water 1/k); c = sqrt(g L_k)
        const double Lk = (seeded && kw>0.0) ? tanh(kw*h)/kw : h;
        const double c = sqrt(g*Lk);
        const double tstar = p->A503*sqrt(Lk/g);
        
        // initial threshold eta_t^I / c: A 554; a new event of a seed (A 550 4) starts at the deta/dt of the column if
        // that is larger, so that B starts from 0 and the viscosity grows continuously
        if(age<=0.0 || r0<=0.0)
        r0 = (seedc && c>0.0) ? MAX(p->A554, et/c) : p->A554;
        
        const double etI = r0*c;
        const double etF = p->A501*c;
        
        const double ets = (age>0.0 && tstar>0.0) ? etI + MIN(age/tstar,1.0)*(etF - etI) : etI;
        
            if(et>ets && ets>0.0)
            {
            const double B = (et>=2.0*ets) ? 1.0 : et/ets - 1.0;
            
            d->vb(i,j) = B*p->A504*p->A504*Lk*et;
            brkage(i,j) = age + p->dt;
            brkk(i,j) = kw;
            brkr(i,j) = r0;
            brkflag(i,j) = 1;
            ++count;
            }
            
            // A 550 4: the event of a seeded column ages also below the threshold, so the threshold relaxes from
            // eta_t^I at the onset of the crest and the viscosity grows continuously over T*
            else
            if(seedc)
            {
            brkage(i,j) = age + p->dt;
            brkk(i,j) = kw;
            brkr(i,j) = r0;
            }
        }
    }
    
    // spread over the band of A 502 h: bd is the width in cells, by marks the band columns
    int nband=0;
    
    SLICELOOP4
    {
    double dxl = p->DXN[IP];
    
    if(p->j_dir==1)
    dxl = MIN(dxl,p->DYN[JP]);
    
    bd(i,j) = (p->wet[IJ]==1 && brkzone(i,j)==0) ? int(ceil(p->A502*MAX(WL(i,j),0.0)/dxl - 1.0e-9)) : 0;
    nband = MAX(nband,bd(i,j));
    }
    
    nband = pgc->globalimax(nband);
    
    SLICELOOP4
    d->vbk(i,j) = brkk(i,j);
    
    for(int qn=0; qn<nband; ++qn)
    {
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4(p,d->vbk,1);
    
    SLICELOOP4
    {
    brkvb(i,j) = d->vb(i,j);
    bt_vbs(i,j) = d->vbk(i,j);
    
        if(qn<bd(i,j))
        {
        auto take = [&](int ii, int jj) { if(d->vb(ii,jj)>brkvb(i,j)) { brkvb(i,j) = d->vb(ii,jj); bt_vbs(i,j) = d->vbk(ii,jj); } };
        take(i-1,j);
        take(i+1,j);
        
        if(p->j_dir==1)
        {
        take(i,j-1);
        take(i,j+1);
        }
        }
    }
    
    SLICELOOP4
    {
    d->vb(i,j) = brkvb(i,j);
    d->vbk(i,j) = bt_vbs(i,j);
    }
    }
    
    SLICELOOP4
    d->breaking(i,j) = brkflag(i,j);
    
    pgc->gcsl_start4(p,d->vb,1);
    pgc->gcsl_start4(p,d->vbk,1);
    pgc->gcsl_start4int(p,d->breaking,50);
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}

// the relaxation zones of the wave generation and the numerical beach (ioflow::relax_boxes): deta/dt of the
// continuity is not the change of the relaxed water level there
void nhflow_breaking::breaking_zones(lexer *p, ioflow *pflow)
{
    if(brkzone_built || p->A550<2 || pflow==nullptr)
    return;
    
    // boxes xs, xe, ys, ye
    std::vector<double> b;
    pflow->relax_boxes(p,b);
    
    SLICELOOP4
    {
    const double x = (p->geo_curv==1) ? p->XC2D[(i-p->imin)*p->jmax + (j-p->jmin)] : p->XP[IP];
    const double y = (p->geo_curv==1) ? p->YC2D[(i-p->imin)*p->jmax + (j-p->jmin)] : p->YP[JP];
    
    brkzone(i,j) = 0;
    
        for(size_t q=0; q+3<b.size(); q+=4)
        if(x>=b[q] && x<=b[q+1] && (p->j_dir==0 || (y>=b[q+2] && y<=b[q+3])))
        brkzone(i,j) = 1;
    }
    
    brkzone_built = true;
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
    // deta/dt of the stage from the continuity (nhflow_fsf rk2/rk3 steps). (eta - eta_n)/(alpha dt) was right
    // only in the first stage: eta_n is the level of the previous stage, so in the RK2 corrector it was about
    // dt d2eta/dt2 and the test did not fire
    const double g = fabs(p->W22);
    
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(d->detadt(i,j) > p->A554*sqrt(g*MAX(WL(i,j),0.0)))
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
        // the cell and its direct neighbours A 557, the next ones A 557/2 (before, the second test
        // overwrote the direct neighbours with A 557/2)
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0)
        d->vb(i,j) = p->A557;
            
        else
        if(brkflag(i-2,j)>0 || brkflag(i+2,j)>0)
        d->vb(i,j) = 0.5*p->A557;
    }

    if(p->j_dir==1)
    SLICELOOP4
    {   
        if(brkflag(i,j)>0 || brkflag(i-1,j)>0 || brkflag(i+1,j)>0 || brkflag(i,j-1)>0 || brkflag(i,j+1)>0)
        d->vb(i,j) = p->A557;
            
        else
        if(( brkflag(i-1,j-1)>0 || brkflag(i-1,j+1)>0 || brkflag(i+1,j-1)>0 || brkflag(i+1,j+1)>0
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
        
    // the flags for the output (VTP 'breaking') and the breaking count of the log
    SLICELOOP4
    d->breaking(i,j) = brkflag(i,j);
    
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
