/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_force.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_force::nhflow_force(lexer* p, fdm_nhf *d, ghostcell *pgc, int qn) : interfac(1.6),zero(0.0),ID(qn)
{
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_NHFLOW_SOLID",0777);
	
    forceprintcount=0;
    
    p->Darray(eta,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(vertice,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(nodeflag,p->imax*p->jmax*(p->kmax+2));
    
    // open files
    print_ini(p,d,pgc);
    
    is = p->posc_i(p->P81_xs[ID]);
    ie = p->posc_i(p->P81_xe[ID]);
    
    js = p->posc_j(p->P81_ys[ID]);
    je = p->posc_j(p->P81_ye[ID]);
    
    ks = p->posc_k(p->P81_zs[ID]);
    ke = p->posc_k(p->P81_ze[ID]);
	
	xs = p->P81_xs[ID];
	xe = p->P81_xe[ID];
	
	ys = p->P81_ys[ID];
	ye = p->P81_ye[ID];
	
	zs = p->P81_zs[ID];
	ze = p->P81_ze[ID];
	
	xm = xs + (xe-xs)*0.5;
	ym = ys + (ye-ys)*0.5;
	zm = zs + (ze-zs)*0.5;
	
    gcval_press=40;  
}

nhflow_force::~nhflow_force()
{
}

void nhflow_force::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    triangulation(p,d,pgc);
	reconstruct(p,d);
	print_vtp(p,d,pgc);
} 

void nhflow_force::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	// forcecalc
    triangulation(p,d,pgc);
	reconstruct(p,d);
    force_calc(p,d,pgc);

        if(p->mpirank==0)
        {
        if(p->count==2)
        cout<<"Atot_solid: "<<A_tot<<endl;  
        
        cout<<"Ax: "<<Ax<<" Ay: "<<Ay<<" Az: "<<Az<<endl;
        cout<<"Atot: "<<A_tot<<endl;
        cout<<"Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<endl;

        print_force(p,d,pgc);
        }
        
    print_vtp(p,d,pgc);
    
    deallocate(p,d,pgc);
} 

