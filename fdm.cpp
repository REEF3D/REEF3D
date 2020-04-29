/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fdm.h"
#include"lexer.h"

fdm::fdm(lexer *p)
			:u(p),F(p),
			v(p),G(p),
			w(p),H(p),
			press(p),
            Fi(p),
			eddyv(p),
			L(p),
			ro(p),visc(p),visctot(p),
			phi(p),vof(p),
			conc(p),
            topo(p),solid(p),
            test(p),
			fb(p),porosity(p),
			walld(p),
			nodeval(p),flag(p),etaloc(p),
            eta(p),eta_n(p),WL(p),depth(p),Fifsf(p),
            Fz(p),
            bed(p),bedzh(p),bedk(p),wet(p),
            bedload(p),qbx(p),qby(p),
            P(p),Q(p),K(p),
			xvec(p),rhsvec(p),M(p)
{

	

	maxF=0.0;
	maxG=0.0; 
	maxH=0.0;
	maxK=0.0;
	maxE=0.0;
    
	gi=p->W20;
	gj=p->W21;
	gk=p->W22;


    p->Iarray(pvccnode,p->facetnum*4,8);
	p->Iarray(ccedge,p->facetnum*4);
    
    
    C1.allocate(p);
    C2.allocate(p);
    C3.allocate(p);
    C4.allocate(p);
    C4a.allocate(p);
    C6.allocate(p);
}


















