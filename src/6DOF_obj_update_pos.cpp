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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>

void sixdof_obj_cfd::update_position_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalize)
{
    // Calculate new position
    rb.euler_angles();

    // Update STL mesh
    update_trimesh_3D(p,a,pgc,finalize);
    
    // Update angular velocities 
    rb.update_omega();
    
    // the pose and velocity: in the 6DOF log files (REEF3D_*_6DOF), no longer on the terminal

}

void sixdof_obj_cfd::update_trimesh_3D(lexer *p, fdm *a, ghostcell *pgc, bool finalize)
{
	// Update position of triangles 
    geom.transform(R_,c_);

    // the triangles as at the last cast on this grid (a fixed body, or a body at rest): the level set
    // of that cast. The ray cast and the reinitialisation of a fixed hull took more than half of the
    // run time at every stage (OF1 in current); the level set is the same to the bit.
    if(fb_unchanged(p,a,pgc))
    {
        std::copy(fbc_v.begin(),fbc_v.end(),a->fb.V);
        return;
    }

    // Update floating level set function
	ray_cast(p,a,pgc);
	reini_RK2(p,a,pgc,a->fb);
    
    pgc->start4a(p,a->fb,50);
    
    fb_store(p,a);
}

// true on all ranks if the triangles of the body are those of the last cast on the grid of a (bitwise),
// and that cast is stored; deforming (X 195) or remeshed bodies (X 185) have other triangles
bool sixdof_obj_cfd::fb_unchanged(lexer *p, fdm *a, ghostcell *pgc)
{
    const int n = amr_tricount();
    int same = (fbc_a==a && (int)fbc_v.size()==p->imax*p->jmax*p->kmax && (int)fbc_tri.size()==9*n) ? 1 : 0;
    
    size_t c = 0;
    for(int d=0; d<3 && same; ++d)
    {
        double **t = amr_tri(d);
        for(int q=0; q<n && same; ++q)
        for(int e=0; e<3; ++e)
        {
            if(fbc_tri[c]!=t[q][e])
            same = 0;
            ++c;
        }
    }
    
    // the ray cast and the reinitialisation exchange ghost cells: all ranks take the same branch
    return pgc->globalimin(same)==1;
}

void sixdof_obj_cfd::fb_store(lexer *p, fdm *a)
{
    const int n = amr_tricount();
    
    fbc_tri.resize(9*n);
    size_t c = 0;
    for(int d=0; d<3; ++d)
    {
        double **t = amr_tri(d);
        for(int q=0; q<n; ++q)
        for(int e=0; e<3; ++e)
        fbc_tri[c++] = t[q][e];
    }
    
    fbc_v.assign(a->fb.V, a->fb.V + p->imax*p->jmax*p->kmax);
    fbc_a = a;
}

