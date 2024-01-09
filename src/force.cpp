/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"force.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

force::force(lexer* p, fdm *a, ghostcell *pgc, int qn):nodefill(p),vertice(p),nodeflag(p),interfac(1.6),zero(0.0),eta(p),ID(qn)
{
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SOLID",0777);
	
    forceprintcount=0;
    
    
    // open files
    print_ini(p,a,pgc);
    
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

force::~force()
{
}
void force::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    triangulation(p,a,pgc,a->phi);
	reconstruct(p,a,a->phi);
    pgc->gcxsd_seed(p,a);
    pgc->gcbsd_seed(p,a);
	
	print_vtp(p,a,pgc);
} 

void force::start(lexer *p, fdm *a, ghostcell *pgc)
{
    pgc->start4(p,a->press,gcval_press);
    pgc->gcxsd_update(p, a, a->press);
    pgc->gcbsd_update(p, a, a->press);
    
	// forcecalc
    force_calc(p,a,pgc);
    
        if(p->mpirank==0)
        {
        if(p->count==2)
        cout<<"Atot_solid: "<<A_tot<<endl;  
        
        cout<<"Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<endl;

        print_force(p,a,pgc);
        }
    
    pgc->start4(p,a->press,gcval_press);
} 

