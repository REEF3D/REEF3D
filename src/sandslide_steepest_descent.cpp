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

#include"sandslide_steepest_descent.h"
#include"sediment_fdm.h"
#include"sediment_mixture.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sliceint.h"

sandslide_steepest_descent::sandslide_steepest_descent(lexer *p) : norm_vec(p), bedslope(p), fh(p)
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
}

sandslide_steepest_descent::~sandslide_steepest_descent()
{
}

void sandslide_steepest_descent::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
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
        
        if(s->pmix!=nullptr)
        s->pmix->slide_zero(p,pgc);
        
        // slide loop
       slide(p,pgc,s);

        
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

        count=pgc->globalimax(count);

        p->slidecells=count;
        

        if(p->mpirank==0)
        cout<<"sandslide_steepest_descent corrections: "<<p->slidecells<<endl;
        
        if(p->slidecells==0)
        break;
    }
}

void sandslide_steepest_descent::slide(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double tan_phi;
    
        // Reset flux accumulator
    SEDSLICELOOP
    fh(i,j) = 0.0;
    
    SEDSLICELOOP
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    {
        int i_steep, j_steep;
        double max_slope, dist_steep;
        
        tan_phi = tan(s->phi(i,j));
        
        // Find the steepest downslope neighbor
        find_steepest_neighbor(p, s, s->bedzh, i, j, i_steep, j_steep, max_slope, dist_steep);
        
        // Check if slope exceeds angle of repose
        if(max_slope > tan_phi)
        {
            double z0 = s->bedzh(i,j);
            double z1 = s->bedzh(i_steep, j_steep);
            double dz = z0 - z1;
            
            // Equilibrium elevation difference for this distance
            double dz_eq = tan_phi * dist_steep;
            
            // Excess height above equilibrium
            double excess = dz - dz_eq;
            
            // Amount to transfer: half the excess moves to achieve equilibrium
            // (because both cells adjust: source lowers, target rises)
            // Apply relaxation factor for stability
            double transfer = relax * 0.5 * excess;
            
            // Accumulate flux (don't modify topo directly yet)
            fh(i,j) -= transfer;
            fh(i_steep, j_steep) += transfer*SLIDE_AR(i_steep-i,j_steep-j);
            
            if(s->pmix!=nullptr)
            s->pmix->slide_transfer(i,j,i_steep,j_steep,transfer);
            
            ++count;
        }
    }
        
}


void sandslide_steepest_descent::find_steepest_neighbor(lexer *p,sediment_fdm *s, slice &zh, int i, int j, int& i_steep, int& j_steep,
                                                        double& max_slope, double& dist_steep)
{
    double z0 = zh(i,j);
    max_slope = 0.0;
    i_steep = i;
    j_steep = j;
    dist_steep = p->DXP[IP];
    
    // 8-connectivity: check all surrounding cells
    for(int di = -1; di <= 1; ++di)
    {
        for(int dj = -1; dj <= 1; ++dj)
        {
            if((di == 0 && dj == 0)||s->DFBED[(i-p->imin+di)*p->jmax + (j-p->jmin+dj)]<0)
            continue;
            
            // no transfer into physical boundary ghost cells (lost) or across rows in 2D
            if(!SLIDE_NB(di,dj))
            continue;
                
            // distance between the cell centres (was 0.5*(dx+dy) and sqrt(2) of it: wrong for dx != dy)
            double ddx = di<0?p->DXP[IM1]:(di>0?p->DXP[IP]:0.0);
            double ddy = dj<0?p->DYP[JM1]:(dj>0?p->DYP[JP]:0.0);
            double dist = sqrt(ddx*ddx + ddy*ddy);
                
            // Elevation difference (positive means neighbor is lower)
            double dz = z0 - zh(i+di, j+dj);
                
            // Slope in this direction
            double slope = dz / dist;
                
            // Track steepest downslope direction
            if(slope > max_slope)
            {
                max_slope = slope;
                i_steep = i + di;
                j_steep = j + dj;
                dist_steep = dist;
            }
        }
    }
}
