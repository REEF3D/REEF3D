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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj_nhflow::update_position_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &fsglobal, bool finalize)
{
    // Calculate new position
    rb.euler_angles();
    
    // Update STL mesh
    update_trimesh_nhflow(p,d,pgc,finalize);

    // Update angular velocities 
    rb.update_omega();
    
    k=p->knoz-1;
    
    SLICELOOP4
    fsglobal(i,j) = d->FB[IJK];
    
    pgc->gcsl_start4(p,fsglobal,50);
    
    if(p->mpirank==0 && finalize==true)
    {
        cout<<"XG: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
        cout<<"Ue: "<<u_fb(0)<<" Ve: "<< u_fb(1)<<" We: "<< u_fb(2)<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }
}

void sixdof_obj_nhflow::update_trimesh_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, bool finalize)
{
	// Update position of triangles 
    geom.transform(R_,c_);
    
    // Update floating level set function
	ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);
}

void sixdof_obj_nhflow::update_wavemaker_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &fsglobal, int iter, bool finalize)
{
    // paddle at the time of the stage result, like the floating body: stage s of TVD RK3 gives
    // the solution at t^n + (1, 1/2, 1) dt, of RK2 at t^n + (1, 1) dt (p->simtime is t^n)
    static const double crk3[3] = {1.0, 0.5, 1.0};
    static const double crk2[2] = {1.0, 1.0};
    
    double t = p->simtime + p->dt;
    
    if(p->A510==3 && iter>=0 && iter<3)
    t = p->simtime + crk3[iter]*p->dt;
    
    if(p->A510==2 && iter>=0 && iter<2)
    t = p->simtime + crk2[iter]*p->dt;
    
    // Update position of wavemaker
    if(p->X170==1)
	piston(p,pgc,0,t);
    
    if(p->X172==1)
    flap_double(p,pgc,0,t);
    
    // Update floating level set function
	ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);

}