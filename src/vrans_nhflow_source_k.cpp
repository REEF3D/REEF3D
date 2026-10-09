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

// Porous-media turbulence closure (Nakayama & Kuwahara 1999), B 295 1:
//   k_inf = 3.7 (1-n) n^1.5 |u|^2,   eps_inf = 39 (1-n)^2.5 n^2 |u|^3 / d50
// u = Darcy velocity, n = static-structure porosity n_s, d50 of the cell (porous_coeff, B 202 layers). Sources are explicit, added to the
// rows of the k / eps / omega systems (row order = LOOP order, one row per cell).
void vrans_nhflow_f::turb_inf(lexer *p, fdm_nhf *d, double &kinf, double &einf)
{
    double uu;
    
    uu = d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK]*p->y_dir + d->W[IJK]*d->W[IJK];
    
    kinf = 3.7*(1.0-porval)*pow(porval,1.5)*uu;
    einf = 39.0*pow(1.0-porval,2.5)*pow(porval,2.0)*pow(uu,1.5)/partval;
    
    kinf = MAX(kinf,0.0);
    einf = MAX(einf,0.0);
}

// k-epsilon, k-equation: + n eps_inf
void vrans_nhflow_f::ke_source(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    double kinf,einf;
    
    count=0;
    if(p->B295==1)
    LOOP
    {
        porous_coeff(p,d);
        
        if(H>1.0e-12)
        {
        turb_inf(p,d,kinf,einf);
        
        d->rhsvec.V[count] += porval*einf;
        }
        
    ++count;
    }
}

// k-omega, k-equation: + n beta* k_inf omega_inf = + n eps_inf   (omega_inf = eps_inf/(cmu k_inf))
void vrans_nhflow_f::kw_source(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    double kinf,einf;
    
    count=0;
    if(p->B295==1)
    LOOP
    {
        porous_coeff(p,d);
        
        if(H>1.0e-12)
        {
        turb_inf(p,d,kinf,einf);
        
        d->rhsvec.V[count] += porval*einf;
        }
        
    ++count;
    }
}
