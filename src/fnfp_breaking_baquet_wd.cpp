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

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_breaking::breaking_baquet_wd(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    int ii,jj;
    
    if(p->count>count_n)
    {
    SLICELOOP4
    c->breaking(i,j)=0;
    
    count_n=p->count;
    }
    
    SLICELOOP4
    {
    bx(i,j)=0;
    by(i,j)=0;
    }
    
    pgc->gcsl_start4int(p,c->breaking,50);
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    
    // steepness induced breaking
    if((p->A351==2 || p->A351==3) && p->count>1)
    {
    SLICELOOP4
    {
            // x
            if(c->Ex(i,j)   < -p->A355)
            bx(i,j) = 10;
            
            if(c->Ex(i,j)   > p->A355)
            bx(i,j) = 20;

            // y
            if(p->j_dir==1)
            if( c->Ey(i,j)   < -p->A355)
            by(i,j) = 10;
            
            if(p->j_dir==1)
            if( c->Ey(i,j)   > p->A355)
            by(i,j) = 20;
    }
    
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    
    
    
    SLICELOOP4
    if(bx(i,j)>0 || by(i,j)>0)
    {
    c->breaking(i,j)=1;
    }
    
    SLICELOOP4
    {
    // x
    if(bx(i,j)==10 || bx(i-1,j)==10 || bx(i+1,j)==10 || bx(i+2,j)==10)
    c->breaking(i,j)=1;
    
    if(bx(i,j)==20 || bx(i-1,j)==20 || bx(i-2,j)==20 || bx(i+1,j)==20)
    c->breaking(i,j)=1;
    
    // y
    if(by(i,j)==10 || by(i,j-1)==10 || by(i,j+1)==10 || by(i,j+2)==10)
    c->breaking(i,j)=1;
    
    if(by(i,j)==20 || by(i,j-1)==20 || by(i,j-2)==20 || by(i,j+1)==20)
    c->breaking(i,j)=1;
        
    }
    
    }
    
    
    
    SLICELOOP4
    bx(i,j)=0;
    
    pgc->gcsl_start4int(p,bx,50);
    
    // depth induced breaking
    if((p->A351==1 || p->A351==3) && p->count>1)
    {
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if((eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A354*sqrt(9.81*c->WL(i,j)))
        bx(i,j)=2;
    }
    
    pgc->gcsl_start4int(p,bx,50);
    
    SLICELOOP4
    {
    // x
    if(bx(i,j)==2 || bx(i-1,j)==2 || bx(i-2,j)==2 || bx(i+1,j)==2 || bx(i+2,j)==2)
    c->breaking(i,j)=2;
    
    // y
    if(p->j_dir==1)
    if(bx(i,j)==2 || bx(i,j-1)==2 || bx(i,j-2)==2 || bx(i,j+1)==2 || bx(i,j+2)==2)
    c->breaking(i,j)=2;    
    }
    }
    
    
    // ------------------------
    // fill breaking viscosity
    // ------------------------
    
    SLICELOOP4
    c->vb(i,j) = 0.0;
        
        // coastline viscosity
        /*if(p->A343==1)
        SLICELOOP4
        {
            
            if(c->coastline(i,j)>=0.0 && p->A346>0.0)
            {
                db = c->coastline(i,j);
                
                if(db<dist3)
                {
                c->vb(i,j) = rb3(p,db)*p->A346;
                
                filter(p,c,pgc,eta,p->A361,p->A362);
                filter(p,c,pgc,Fifsf,p->A361,p->A362);
                }
            }
        }*/
        
    if(p->j_dir==0)
    SLICELOOP4
    {   
        if(c->breaking(i,j)>0 || c->breaking(i-1,j)>0 || c->breaking(i+1,j)>0)
        c->vb(i,j) = p->A365;
            
        if(c->breaking(i,j)==0 && 
        (c->breaking(i-1,j)>0 || c->breaking(i+1,j)>0 || c->breaking(i-2,j)>0 || c->breaking(i+2,j)>0))
        c->vb(i,j) = 0.5*p->A365;
    }

    if(p->j_dir==1)
    SLICELOOP4
    {   
        if(c->breaking(i,j)>0 || c->breaking(i-1,j)>0 || c->breaking(i+1,j)>0 || c->breaking(i,j-1)>0 || c->breaking(i,j+1)>0)
        c->vb(i,j) = p->A365;
            
        if(c->breaking(i,j)==0 && 
        ( c->breaking(i-1,j-1)>0 || c->breaking(i-1,j+1)>0 || c->breaking(i+1,j-1)>0 || c->breaking(i+1,j+1)>0
        || c->breaking(i-2,j)>0 || c->breaking(i+2,j)>0 || c->breaking(i,j-2)>0 || c->breaking(i,j+2)>0))
        c->vb(i,j) = 0.5*p->A365;
    }
        
        // coastline / wet-dry front viscosity A346 (A343 1 with A337, A343 2/3 with A332)
        coast_viscosity(p,c,pgc);
        
        // A352 5: the zone viscosity in the cells of the depth-induced criterion instead of the filters below
        if(p->A352==5)
        zone_viscosity(p,c,pgc);
        
        // additional breaking filter
        // shallow
        if(p->A352==1)
        SLICELOOP4
        if(c->breaking(i,j)==2)
        {
         filter(p,c,pgc,eta,p->A361,p->A362);
         filter(p,c,pgc,Fifsf,p->A361,p->A362);
        }   
        
        // deep
        if(p->A352==2)
        SLICELOOP4
        if(c->breaking(i,j)==1)
        {
         filter(p,c,pgc,eta,p->A361,p->A362);
         filter(p,c,pgc,Fifsf,p->A361,p->A362);
        }   
        
        // all
        if(p->A352==3)
        SLICELOOP4
        if(c->breaking(i,j)>=1)
        {
         filter(p,c,pgc,eta,p->A361,p->A362);
         filter(p,c,pgc,Fifsf,p->A361,p->A362);
        }   
        
        
    pgc->gcsl_start4(p,c->vb,1);
    
    
    SLICELOOP4
    c->breaklog(i,j)=0;
    
    // breaklog
    int count=0;
    
    SLICELOOP4
    if(c->breaking(i,j)>0)
    {
    c->breaklog(i,j)=1;
    ++count;
    }

    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}

// coastline / wet-dry front viscosity (raises vb to at least A346 there) with A358 1 (default 0: off, as before
// A337); part of the Baquet model (A350 1), also used by A350 4/5
void fnpf_breaking::coast_viscosity(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // A343 2/3: viscosity band along the moving wet-dry front
    // (replaces the static coastline viscosity above, which is tied to
    // the initial coastline); cells within A332 cells of a dry cell get
    // at least A346, applied through the implicit damping solve
    if(p->A358==1 && p->A343>=2 && p->A332>0 && p->A346>0.0)
    {
        SLICELOOP4
        wdband(i,j) = (p->wet[IJ]==0) ? 1 : 0;
        
        pgc->gcsl_start4int(p,wdband,50);
        
        for(int qn=0; qn<p->A332; ++qn)
        {
            // tag qn+2 marks cells added in this sweep, so the band grows
            // by exactly one cell per sweep
            SLICELOOP4
            if(wdband(i,j)==0)
            {
                int hit=0;
                
                for(int di=-1; di<=1; ++di)
                for(int dj=-p->j_dir; dj<=p->j_dir; ++dj)
                {
                    const int b = wdband(i+di,j+dj);
                    
                    if(b>0 && b<=qn+1)
                    hit=1;
                }
                
                if(hit==1)
                wdband(i,j) = qn+2;
            }
            
            pgc->gcsl_start4int(p,wdband,50);
        }
        
        SLICELOOP4
        if(p->wet[IJ]==1 && wdband(i,j)>1)
        c->vb(i,j) = MAX(c->vb(i,j), p->A346);
    }
    
    // A343 1: coastline viscosity A346 through the implicit damping solve (needs A350>0), relaxed from zero at
    // the edge of the coastline band (A341 cells or A342 m, dist3) to A346 at the coastline with the weight
    // 1 - rb3, the complement of the relaxation of eta and Fi in the band: no viscosity without a band (A341 0,
    // A342 0, the default) or with A346 0, no on/off switch at the band edge
    //   A337 0: legacy, no coastline viscosity (as master, where the block was commented out)
    //   A337 1: the ramp in the band
    //   A337 2, 3: the ramp, and A346 in the thin-film cells (WL < 10*A344), where a drawdown has exposed the
    //           bed but the static mask keeps the cell wet and clamped
    if(p->A358==1 && p->A343==1 && p->A337>=1 && p->A346>0.0)
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(c->coastline(i,j)>=0.0 && dist3>0.0 && c->coastline(i,j)<dist3)
        c->vb(i,j) = MAX(c->vb(i,j), (1.0-rb3(p,c->coastline(i,j)))*p->A346);

        if(p->A337>=2 && c->WL(i,j) < 10.0*c->wd_criterion)
        c->vb(i,j) = MAX(c->vb(i,j), p->A346);
    }
}

// A352 5 (A350 1, 4, 5): instead of the filter the breaking viscosity is raised to at least nu = dx c/(4 C), C = 0.5,
// in the cells of the depth-induced criterion (flag 2: eta_t > A354 sqrt(g h), +-2 cells) and +-2 cells around them,
// c the phase speed of the 2 dx wave in the local depth: a 2 dx spike decays by e while it travels half a cell,
// longer waves much less (about with the square of their length). The maximum with the breaking and coastline
// viscosity, applied by fnpf_fsfbc::damping. With A350 4/5 the flag-2 cells are not limited to the still-water
// depth of the deepest onset (bt_db): at a steep shore the crest tracking sees the onset late and shallow
// (FNPF Breaking validation 09).
void fnpf_breaking::zone_viscosity(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const double g = fabs(p->W22);
    const double C = 0.5;

    SLICELOOP4
    bx(i,j) = (c->breaking(i,j)>=2 && p->wet[IJ]==1) ? 1 : 0;

    pgc->gcsl_start4int(p,bx,50);

    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        int hit=0;

        for(int qi=-2; qi<=2; ++qi)
        for(int qj=-2*p->j_dir; qj<=2*p->j_dir; ++qj)
        if(bx(i+qi,j+qj)==1)
        hit=1;

        if(hit==1)
        {
            const double wl = MAX(c->WL(i,j), c->wd_criterion);
            const double dx = (p->j_dir==1) ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];
            const double k = PI/dx;
            const double cp = sqrt(g*tanh(k*wl)/k);

            c->vb(i,j) = MAX(c->vb(i,j), dx*cp/(4.0*C));
        }
    }

    pgc->gcsl_start4(p,c->vb,1);
}
