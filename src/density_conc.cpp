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

#include"interface_width.h"
#include"density_conc.h"
#include"lexer.h"
#include"fdm.h"
#include"concentration.h"
#include"heaviside.h"

density_conc::density_conc(lexer* p, concentration *& ppconc) 
{
        pconc = ppconc;
    
        
        H=0.0;
}

density_conc::~density_conc()
{
}

double density_conc::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double concval;

        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
        
        concval = 0.5*(pconc->val(i,j,k) + pconc->val(i+aa,j+bb,k+cc));
        

        H = heaviside(phival,interface_width_face(p,a->phi,i,j,k,aa,bb,cc));
        
        roval = (p->W1+concval*p->C1)*H + (p->W3+concval*p->C3)*(1.0-H);
    

	
	return roval;		
}




