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


#include"6DOF_obj.h"
#include"6DOF_obj_fnpf.h"
#include"lexer.h"
#include"ghostcell.h"

//  Mesh refinement with subcycling (G 7 1: nhflow_amr, fnpf_amr): the finest level advances the
//  body; the coarser levels step before it with a predicted copy of it.  amr_save keeps the
//  rigid-body state, amr_restore puts it back with the trimesh and the body velocities.

void sixdof_obj::amr_save()
{
    rb_amr = rb;
}

void sixdof_obj::amr_restore(lexer *p, ghostcell *pgc)
{
    rb = rb_amr;
    geom.transform(R_,c_);
    update_fbvel(p,pgc);
}

// FNPF: the RK stage of the predicted copy with the loads and the added mass of the last solve
void sixdof_obj_fnpf::amr_predict_fnpf(lexer *p, ghostcell *pgc, int iter)
{
    if(p->A310==3)
    rk3(p,pgc,iter);
    
    else
    rk4(p,pgc,iter);
    
    update_position_fnpf(p,pgc,false);
}
