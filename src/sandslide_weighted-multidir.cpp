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

#include"sandslide_weighted_multidir.h"
#include"sediment_fdm.h"
#include"sediment_mixture.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sliceint.h"

sandslide_weighted_multidir::sandslide_weighted_multidir(lexer *p) : norm_vec(p), bedslope(p), fh(p)
{
    if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	fac1 = p->S92*(1.0/6.0);
	fac2 = p->S92*(1.0/12.0);
    
    // S 92: correction factor of the transfers (1: each over-steep pair alone is brought exactly to phi)
    relax=p->S92;
    
    // converged when no neighbour pair is steeper than phi by more than tol (height)
    tol = 1.0e-4*MAX(p->S20,1.0e-6);
    
}

sandslide_weighted_multidir::~sandslide_weighted_multidir()
{
}

void sandslide_weighted_multidir::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    SEDSLICELOOP
    s->slide_fh(i,j)=0.0;
    
    // mainloop
    for(int qn=0; qn<p->S91; ++qn)
    {
        slidecount=0;
        
        // fill
        SEDSLICELOOP
        fh(i,j)=0.0;
        
        pgc->gcsl_start4(p,fh,1);
        
        if(s->pmix!=nullptr)
        s->pmix->slide_zero(p,pgc);
        
        // slide loop
       compute_fh(p,pgc,s);

        
        pgc->gcslparax_fh(p,fh,4);
        
        // fill back
        SEDSLICELOOP
        {
        s->slide_fh(i,j)+=fh(i,j);
        s->bedzh(i,j)+=fh(i,j);
        }
        

        pgc->gcsl_start4(p,s->bedzh,1);
        
        // multi-fraction bed: sorting of the slid material
        if(s->pmix!=nullptr)
        s->pmix->slide_finish(p,pgc,s);

        slidecount=pgc->globalimax(slidecount);

        p->slidecells=slidecount;
        

        if(p->mpirank==0)
        cout<<"sandslide_weighted_multidir corrections: "<<p->slidecells<<endl;
        
        if(p->slidecells==0)
        break;
    }
}

void sandslide_weighted_multidir::compute_fh(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    // Every cell sends to each lower neighbour k (8-connectivity) whose slope exceeds the angle of
    // repose a share of half the excess height e_k = dz_k - d_k*tan(phi):
    //      t_k = relax * 0.5*e_k * e_k/sum(e)
    // A single over-steep pair gets exactly half its excess (pair equilibrium for equal cells);
    // the total outflow is at most half the largest excess, and every pair keeps draining until
    // it reaches phi.
    // (was: total 0.5*relax*min(excess) with relax 0.5, shared out by excess-slope^2 weights:
    // min(excess) -> 0 as soon as one neighbour is close to phi, so the cell stopped draining its
    // steeper neighbours; the cone of validation case 03 stalled at 39 deg on the diagonals)
    int di_arr[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int dj_arr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    
    auto open = [&](int di, int dj)
    {
        return SLIDE_NB(di,dj) && s->DFBED[(i-p->imin+di)*p->jmax + (j-p->jmin+dj)]>0;
    };
    
    SEDSLICELOOP
    {
        double z0 = s->bedzh(i,j);
        
        double excess[8];
        int ni[8], nj[8];
        double total_excess = 0.0;
        int count = 0;
        int over = 0;
        
        tan_phi = tan(s->phi(i,j));
        
        for(int k = 0; k < 8; ++k)
        {
            int di = di_arr[k];
            int dj = dj_arr[k];
            
            // no transfer into structures, physical boundary ghost cells or across rows in 2D;
            // a diagonal needs at least one of the two cells beside it open
            if(!open(di,dj))
            continue;
            
            if(di!=0 && dj!=0 && !open(di,0) && !open(0,dj))
            continue;
            
            // distance between the cell centres
            double ddx = di<0?p->DXP[IM1]:(di>0?p->DXP[IP]:0.0);
            double ddy = dj<0?p->DYP[JM1]:(dj>0?p->DYP[JP]:0.0);
            double d = sqrt(ddx*ddx + ddy*ddy);
            
            // excess height over the angle of repose (positive: steeper than phi, downslope)
            double e = (z0 - s->bedzh(i+di, j+dj)) - tan_phi*d;
            
            if(e > 0.0)
            {
                ni[count] = i + di;
                nj[count] = j + dj;
                excess[count] = e;
                total_excess += e;
                ++count;
                
                if(e > tol)
                over = 1;
            }
        }
        
        slidecount += over;
        
        for(int k = 0; k < count; ++k)
        {
            double t = relax*0.5*excess[k]*excess[k]/total_excess;
            
            fh(i,j) -= t;
            fh(ni[k], nj[k]) += t*SLIDE_AR(ni[k]-i,nj[k]-j);
            
            if(s->pmix!=nullptr)
            s->pmix->slide_transfer(i,j,ni[k],nj[k],t);
        }
    }
}


double  sandslide_weighted_multidir::compute_slope(lexer* p, slice & zh, int i, int j, int di, int dj)
{
    // Compute distance to neighbor (accounting for diagonal)
    double dist;
    double dx = 0.5*(p->DXN[IP] + p->DYN[JP]);
    
    if(di != 0 && dj != 0)
        dist = dx * sqrt(2.0);  // diagonal
    else
        dist = dx;              // cardinal direction
    
    // Slope = elevation difference / horizontal distance
    double dz = zh(i,j) - zh(i+di, j+dj);
    
    return dz / dist;
}
