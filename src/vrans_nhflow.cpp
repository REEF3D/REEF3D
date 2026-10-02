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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

vrans_nhflow_f::vrans_nhflow_f(lexer *p, fdm_nhf *d, ghostcell *pgc) : nhflow_geometry(p,d,pgc), Cval(p->B264)
{
    p->Darray(UN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WN,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(P,p->imax*p->jmax*(p->kmax+2));
    
    print_force_ini(p,d,pgc);
    
    cmfac=1.0;
}

vrans_nhflow_f::~vrans_nhflow_f()
{
}

void vrans_nhflow_f::update(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, int val)
{
    ray_cast(p, d, pgc, d->PORSTRUC);
    reini_RK2(p, d, pgc, d->PORSTRUC);

    LOOP
    {
    H = Hporface(p,d,0,0,0);   
    
    d->POR[IJK]     = H*p->B201_n + (1.0-H)*1.0;
	d->PORPART[IJK] = H*p->B201_d50;
    }
    
    // porous floating body (X 16): POR is reset above, so re-apply the moving body porosity
    // n = 1 - H_fb(1 - n_fb), same Heaviside as sixdof_obj::Hsolidface_nhflow()
    if(p->X10>0 && p->X16==1)
    LOOP
    {
        double psi,Hfb;
        
        if(p->j_dir==0)
        psi = p->A526*p->DXN[IP];
        
        if(p->j_dir==1)
        psi = p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);
        
        Hfb = 0.5*(1.0 + (-d->FB[IJK])/psi + (1.0/PI)*sin((PI*(-d->FB[IJK]))/psi));
        
        if(-d->FB[IJK] > psi)
        Hfb = 1.0;
        
        if(-d->FB[IJK] < -psi)
        Hfb = 0.0;
        
        d->POR[IJK] = MIN(d->POR[IJK], 1.0 - Hfb*(1.0 - p->X16_n));
    }
    
    // porosity of the DEM particles (E 28), 1 otherwise
    LOOP
    d->POR[IJK] *= d->PORDEM[IJK];

    pgc->start5Vfull(p,d->POR,1);
    pgc->start5Vfull(p,d->PORPART,1);
    
    // print force
    if(p->B208==1)
    force_calc(p, d, pgc, alpha, val);
}

