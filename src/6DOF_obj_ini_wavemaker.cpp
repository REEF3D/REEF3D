/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"momentum.h"
#include"ghostcell.h"
#include<sys/stat.h>

void sixdof_obj::initialize_wavemaker(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &WL)
{
    if(p->mpirank==0)
    cout<<"6DOF_obj_ini_wavemaker "<<p->X170<<endl;
    
    if(p->X170==1)
    read_format_1(p,pgc);
    
    // Initialise folder structure
    if(p->X50==1)
	print_ini_vtp(p,pgc);
    
    if(p->X50==2)
    print_ini_stl(p,pgc);
    
    // Initialise processor boundaries
    ini_parallel(p,pgc);
    
    // Initialise objects
	objects_create(p,pgc);
    
    // Recalculate distances
	ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);
    
    // Print initial body 
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    
}
