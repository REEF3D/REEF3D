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

// X 13 2 (nhflow_fsf_hull): the hull as the top of the water columns.
//
// keel_nhflow: vertical ray through the centre of every column (geo_raycast::column_hits, the
// trimesh at the current pose). A closed body gives an even number of crossings; the lowest one is
// the keel S. The column is covered when it has crossings and the keel lies below the level ZETA
// of the column (the free surface of a free column, the piezometric head of a footprint column). The body velocity at the keel point follows from the rigid-body velocity u_fb
// (translation and rotation about the centre of gravity), as in target_velocity.

void sixdof_obj_nhflow::keel_nhflow(lexer *p, fdm_nhf *d, slice &ZETA, slice &S, slice &cover, slice &us, slice &vs, slice &ws)
{
    SLICELOOP4
    {
    cover(i,j) = 0.0;
    us(i,j) = vs(i,j) = ws(i,j) = 0.0;
    }

    if(tricount<=0)
    return;

    vector<vector<double> > hits;
    georay.column_hits(p,tri_x,tri_y,tri_z,0,tricount,p->DXM,hits);

    const int ny = p->knoy;

    // 2D: the hull is extruded in y, the keel point lies in its mid-plane
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

    ILOOP
    JLOOP
    {
        if(p->wet[IJ]==0)
        continue;

        const vector<double> &h = hits[(size_t)i*ny+j];

        if(h.size()<2)
        continue;

        double zk = 1.0e20;

        for(double z : h)
        zk = MIN(zk,z);

        // keel above the level of the column (the free surface, or the piezometric head of a
        // footprint column): not covered
        if(zk - d->bed(i,j) - d->depth(i,j) >= ZETA(i,j) - 1.0e-10)
        continue;

        const double x = p->XP[IP];
        const double y = (p->j_dir==0) ? ymid : p->YP[JP];

        S(i,j) = zk;
        cover(i,j) = 1.0;

        us(i,j) = u_fb(0) + u_fb(4)*(zk - c_(2)) - u_fb(5)*(y - c_(1));
        vs(i,j) = u_fb(1) + u_fb(5)*(x - c_(0)) - u_fb(3)*(zk - c_(2));
        ws(i,j) = u_fb(2) + u_fb(3)*(y - c_(1)) - u_fb(4)*(x - c_(0));
    }
}

// the loads of the footprint columns (nhflow_fsf_hull::forces, already summed over the ranks) into
// the body, as the tail of force_calc_stl: gravity added, the force file written at the end of the step
void sixdof_obj_nhflow::hull_loads_nhflow(lexer *p, ghostcell *pgc, bool finalize, double Fx, double Fy, double Fz,
                                          double Mx, double My, double Mz)
{
    curr_time = finalize ? p->simtime + p->dt : p->simtime;

    Xe = Fx;
    Ye = Fy;
    Ze = Fz;
    Ke = Mx;
    Me = My;
    Ne = Mz;

    // 2D: only surge, heave and pitch
    if(p->j_dir==0)
    Ye = Ke = Ne = 0.0;

    const double Xe_p = Xe, Ye_p = Ye, Ze_p = Ze;

	// Add gravity force
	Xe += p->W20*Mass_fb;
	Ye += p->W21*Mass_fb;
	Ze += p->W22*Mass_fb;

    if(p->mpirank==0 && finalize==1)
    {
        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<0.0<<" \t "<<0.0<<" \t "<<0.0<<endl;
    }
}
