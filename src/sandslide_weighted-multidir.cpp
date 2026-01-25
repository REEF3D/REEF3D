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
    
    relax=0.5;
    
    // S94: weighting method
    // 1 = linear (proportional to excess)
    // 2 = quadratic (smoother)
    // 3 = slope-based
    weight_method = 2;
    
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
        count=0;
        
        // fill
        SEDSLICELOOP
        fh(i,j)=0.0;
        
        pgc->gcsl_start4(p,fh,1);
        
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

        count=pgc->globalimax(count);

        p->slidecells=count;
        

        if(p->mpirank==0)
        cout<<"sandslide_weighted_multidir corrections: "<<p->slidecells<<endl;
        
        if(p->slidecells==0)
        break;
    }
}

void sandslide_weighted_multidir::compute_fh(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    // Neighbor offsets: 8-connectivity
    int di_arr[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int dj_arr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    
    SEDSLICELOOP
    {
        double z0 = s->bedzh(i,j);
        
        // Storage for over-steep neighbors
        double weights[8];
        double excess[8];
        double dist[8];
        int ni[8], nj[8];
        double total_weight = 0.0;
        double dx = 0.5*(p->DXN[IP] + p->DYN[JP]);
        
        int count = 0;
        
        // Check all 8 neighbors
        for(int k = 0; k < 8; ++k)
        {
            int di = di_arr[k];
            int dj = dj_arr[k];
            
            tan_phi = tan(s->phi(i,j));
            
            // Compute distance
            double d;
            if(di != 0 && dj != 0)
                d = dx * sqrt(2.0);
            else
                d = dx;
            
            // Elevation difference (positive = downslope)
            double dz = z0 - s->bedzh(i+di, j+dj);
            
            // Current slope
            double slope = dz / d;
            
            // Check if slope exceeds angle of repose
            if(slope > tan_phi)
            {
                // Excess height that needs to be redistributed
                // This is the height above what would give exactly tan_phi slope
                double excess_height = dz - tan_phi * d;
                
                // Excess slope for weighting
                double excess_slope = slope - tan_phi;
                
                // Store neighbor info
                ni[count] = i + di;
                nj[count] = j + dj;
                dist[count] = d;
                excess[count] = excess_height;
                weights[count] = compute_weight(excess_slope, weight_method);
                total_weight += weights[count];
                count++;
            }
        }
        
        // Distribute flux to over-steep neighbors
        if(count > 0 && total_weight > 0.0)
        {
            // Total volume to move (use minimum excess to ensure stability)
            double min_excess = excess[0];
            for(int k = 1; k < count; ++k)
                min_excess = MIN(min_excess, excess[k]);
            
            // Apply relaxation factor for stability
            // Factor 0.5 ensures we don't overshoot equilibrium
            double total_flux = relax * 0.5 * min_excess;
            
            // Alternative: use weighted average of excess
            // double total_flux = 0.0;
            // for(int k = 0; k < count; ++k)
            //     total_flux += excess[k] * weights[k] / total_weight;
            // total_flux *= relax * 0.5;
            
            // Flux out of current cell
            fh(i,j) -= total_flux;
            
            // Distribute to neighbors proportionally
            for(int k = 0; k < count; ++k)
            {
                double fraction = weights[k] / total_weight;
                fh(ni[k], nj[k]) += total_flux * fraction;
            }
        }
    }
}


double sandslide_weighted_multidir::compute_weight(double excess_slope, int method)
{
    // Different weighting schemes for flux distribution
    
    switch(method)
    {
        case 1:  // Linear: weight proportional to excess slope
            return excess_slope;
            
        case 2:  // Quadratic: smoother distribution, less channeling
            return excess_slope * excess_slope;
            
        case 3:  // Power law: more aggressive toward steepest
            return pow(excess_slope, 1.5);
            
        default:
            return excess_slope;
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
