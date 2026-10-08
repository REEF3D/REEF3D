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
