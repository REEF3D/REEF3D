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

// k-omega, omega-equation: + n beta omega_inf^2,  omega_inf = eps_inf/(cmu k_inf)
// (balances the beta omega^2 sink at omega = omega_inf)
void vrans_nhflow_f::omega_source(lexer *p, fdm_nhf *d, double *KIN, double *EPS, double beta)
{
    double kinf,einf,winf;
    
    count=0;
    if(p->B295==1)
    LOOP
    {
        porous_coeff(p,d);
        
        if(H>1.0e-12)
        {
        turb_inf(p,d,kinf,einf);
        
        if(kinf>1.0e-20)
        {
        winf = einf/(p->cmu*kinf);
        
        d->rhsvec.V[count] += porval*beta*winf*winf;
        }
        }
        
    ++count;
    }
}
