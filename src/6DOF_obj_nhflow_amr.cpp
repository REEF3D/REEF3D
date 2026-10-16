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
Architect: Hans Bihs
--------------------------------------------------------------------*/


#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<vector>

//  Floating bodies on a subcycled NHFLOW mesh refinement (nhflow_amr, G 7 1, X 10 1/2).
//
//  The finest level advances the body: its RK stages with the loads from the finest grid after
//  each stage, as on a uniform fine grid.  The coarser levels step before it (level 0 first,
//  every level before its children) and take their forcing from a predicted copy: the state is
//  saved, advanced over their stages with the loads frozen at their last value (the RK stages
//  only: no external forces, no output) and put back at the end of their step.  Only the cells
//  under the finer levels see the predicted body; the restriction overwrites them.

// stage iter of the body (p->dt, p->simtime: the step of the level), without the level-0 ray
// cast: the grids cast the hull themselves; predict: the RK stage with the frozen loads
void sixdof_obj_nhflow::amr_stage(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, bool predict)
{
    if(predict || p->X10==2)
    solve_eqmotion_oneway_nhflow(p,pgc,iter,false);
    
    else
    solve_eqmotion_nhflow(p,d,pgc,iter,false);
    
    quat_matrices(p);
    rb.euler_angles();
    geom.transform(R_,c_);
    rb.update_omega();
    update_fbvel(p,pgc);
}

// X 13 (nhflow_fsf_body): signed distance of the top cell centre of every column to the hull at
// the current pose; local to the grid of p (no communication): column parity from the vertical
// ray (geo_raycast::column_hits), distance from the triangles within the band of the column
void sixdof_obj_nhflow::top_phi_nhflow(lexer *p, fdm_nhf *d, slice &WL, double band, slice &phi)
{
    SLICELOOP4
    phi(i,j) = band;
    
    if(tricount<=0)
    return;
    
    vector<vector<double> > hits;
    georay.column_hits(p,tri_x,tri_y,tri_z,0,tricount,p->DXM,hits);
    
    const int ny = p->knoy;
    
    // 2D: the hull is extruded in y; the end caps are not walls of the flow
    double ymid=0.0;
    
    if(p->j_dir==0)
    {
        double ylo=1.0e20, yhi=-1.0e20;
        
        for(int n=0; n<tricount; ++n)
        for(int q=0; q<3; ++q)
        {
        ylo = MIN(ylo,tri_y[n][q]);
        yhi = MAX(yhi,tri_y[n][q]);
        }
        
        ymid = 0.5*(ylo+yhi);
    }
    
    // triangles within the band of each column
    vector<vector<int> > cand((size_t)p->knox*ny);
    
    for(int n=0; n<tricount; ++n)
    {
        if(p->j_dir==0)
        {
            const double ax = tri_x[n][1]-tri_x[n][0], ay = tri_y[n][1]-tri_y[n][0], az = tri_z[n][1]-tri_z[n][0];
            const double bx = tri_x[n][2]-tri_x[n][0], by = tri_y[n][2]-tri_y[n][0], bz = tri_z[n][2]-tri_z[n][0];
            const double cx = ay*bz-az*by, cy = az*bx-ax*bz, cz = ax*by-ay*bx;
            const double cn = sqrt(cx*cx+cy*cy+cz*cz);
            
            if(cn<1.0e-30 || fabs(cy)>0.9*cn)
            continue;
        }
        
        const double xmin = MIN(tri_x[n][0],MIN(tri_x[n][1],tri_x[n][2])) - band;
        const double xmax = MAX(tri_x[n][0],MAX(tri_x[n][1],tri_x[n][2])) + band;
        const double ymin = MIN(tri_y[n][0],MIN(tri_y[n][1],tri_y[n][2])) - band;
        const double ymax = MAX(tri_y[n][0],MAX(tri_y[n][1],tri_y[n][2])) + band;
        
        ILOOP
        if(p->XP[IP]>=xmin && p->XP[IP]<=xmax)
        JLOOP
        if(p->j_dir==0 || (p->YP[JP]>=ymin && p->YP[JP]<=ymax))
        cand[(size_t)i*ny+j].push_back(n);
    }
    
    k = p->knoz-1;
    
    ILOOP
    JLOOP
    {
        if(p->wet[IJ]==0)
        continue;
        
        const double x = p->XP[IP];
        const double y = (p->j_dir==0) ? ymid : p->YP[JP];
        const double z = p->ZP[KP]*WL(i,j) + d->bed(i,j);
        
        double d2 = band*band;
        
        for(int n : cand[(size_t)i*ny+j])
        d2 = MIN(d2, geo_raycast::dist2_tri(x,y,z, tri_x[n][0],tri_y[n][0],tri_z[n][0],
                                                  tri_x[n][1],tri_y[n][1],tri_z[n][1],
                                                  tri_x[n][2],tri_y[n][2],tri_z[n][2]));
        
        // inside: an odd number of crossings above the point
        int above=0;
        
        for(double h : hits[(size_t)i*ny+j])
        if(h>z)
        ++above;
        
        const double dist = MIN(sqrt(d2),band);
        
        phi(i,j) = (above%2==1) ? -dist : dist;
    }
}
