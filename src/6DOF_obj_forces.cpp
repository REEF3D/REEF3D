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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"

void sixdof_obj::hydrodynamic_forces_cfd(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter, bool finalize)
{
    if(p->X60==1)
    forces_stl(p,a,pgc,uvel,vvel,wvel,iter,finalize);
    
    if(p->X60==2)
    forces_lsm(p,a,pgc,uvel,vvel,wvel,iter,finalize);
}

void sixdof_obj::update_forces(lexer *p)
{
    // Forces in inertial system: external loads, linear damping, DOF modes
    const double Fext[6] = {Xext + Xe, Yext + Ye, Zext + Ze, Kext + Ke, Mext + Me, Next + Ne};
    
    rb.assemble_loads(Fext);
    
    if(Ffb_(0)!=Ffb_(0))
    cout<<"Ffb_(0)....###"<<endl;
    
    if(Ffb_(1)!=Ffb_(1))
    cout<<"Ffb_(1)....###"<<endl;
    
    if(Ffb_(2)!=Ffb_(2))
    cout<<"Ffb_(2)....###"<<endl;
    
    
    if(Mfb_(0)!=Mfb_(0))
    cout<<"Mfb_(0)....###"<<endl;
    
    if(Mfb_(1)!=Mfb_(1))
    cout<<"Mfb_(1)....###"<<endl;
    
    if(Mfb_(2)!=Mfb_(2))
    cout<<"Mfb_(2)....###"<<endl;
    
    // FNPF: instantaneous added mass (and implicit PTO terms, X 500 2) on the left-hand side
    if(am_on_ || (pto_on_ && pto_implicit_))
    apply_added_mass(p);
}
