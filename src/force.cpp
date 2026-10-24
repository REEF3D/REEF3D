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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"force.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"amr_grid_lookup.h"
#include<sys/stat.h>
#include<sys/types.h>

force::force(lexer* p, fdm *a, ghostcell *pgc, int qn, bool out):vertice(p),nodeflag(p),interfac(1.6),zero(0.0),eta(p),ID(qn)
{
	// Create Folder
	if(p->mpirank==0 && out)
	mkdir("./REEF3D_SOLID",0777);
	
    forceprintcount=0;
    
    tri = facet = nullptr;
    pt = ccpt = nullptr;
    ls = nullptr;
    confac = numfac = numpt = nullptr;
    numtri_mem = numvert_mem = 0;
    
    
    // open files
    if(out)
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
    free_arrays();
}

// the arrays of the triangulation (allocated where it has pieces)
void force::free_arrays()
{
    auto rows_i = [](int **&f, int n) { if(f!=nullptr) { for(int q=0; q<n; ++q) delete [] f[q]; delete [] f; f = nullptr; } };
    auto rows_d = [](double **&f, int n) { if(f!=nullptr) { for(int q=0; q<n; ++q) delete [] f[q]; delete [] f; f = nullptr; } };
    
    rows_i(tri,numtri_mem);
    rows_d(pt,numvert_mem);
    rows_i(facet,numtri_mem);
    rows_d(ccpt,numtri_mem*4);
    delete [] ls;     ls = nullptr;
    delete [] confac; confac = nullptr;
    delete [] numfac; numfac = nullptr;
    delete [] numpt;  numpt = nullptr;
}

// CFD mesh refinement: level 0
void force::amr_level0(const amr_grid_lookup *pamr)
{
    amr = pamr;
    amr_lev = 0;
}

// CFD mesh refinement: a patch of level lev, the box of its interior (xs, xe, ys, ye, zs, ze)
void force::amr_patch(const amr_grid_lookup *pamr, int lev, const double *box)
{
    amr = pamr;
    amr_lev = lev;
    for(int q=0; q<6; ++q)
    amr_box[q] = box[q];
}

// the piece of the surface with its centre at (x,y,z) belongs to this grid: the finest grid there
bool force::amr_own(double x, double y, double z) const
{
    const int L = amr->leaf_level(x,y,z);
    
    if(amr_lev==0)
    return L==0;
    
    return L==amr_lev && x>=amr_box[0] && x<amr_box[1] && y>=amr_box[2] && y<amr_box[3] && z>=amr_box[4] && z<amr_box[5];
}

// a patch: its pieces, added to F (Fx, Fy, Fz, area, Ax, Ay, Px), rank-local
void force::patch_sum(lexer *p, fdm *a, double *F)
{
    force_sum(p,a,F);
}

void force::ini(lexer *p, fdm *a, ghostcell *pgc, bool printvtp)
{
    triangulation(p,a,pgc,a->phi);
	reconstruct(p,a,a->phi);
	
    if(printvtp)
	print_vtp(p,a,pgc);
    else
    {
    polygon_num=facount;
    polygon_sum=0;
    for(n=0;n<polygon_num;++n)
    polygon_sum+=numpt[n];
    vertice_num = ccptcount;
    }
} 

void force::start(lexer *p, fdm *a, ghostcell *pgc)
{
    pgc->start4(p,a->press,gcval_press);

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

