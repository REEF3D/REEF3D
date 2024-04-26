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

#include"fdm2D.h"
#include"lexer.h"

fdm2D::fdm2D(lexer *p)
			:eta(p),eta_n(p),
            P(p),Pn(p),Q(p),Qn(p),
            F(p),G(p),L(p),
            ws(p),
            press(p),
            eddyv(p),kin(p),eps(p),
            bed(p),bed0(p),depth(p),
            solidbed(p),topobed(p),
            bednode(p),
			 hx(p),hy(p),hp(p),
			 xvec(p),rhsvec(p),M(p),
            dpx(p),dpy(p),test(p),Hs(p),fs(p),
            breaking(p),breaking_print(p),
            wet1(p),wet2(p),
			 nodeval(p),
			 cmu(0.09),
            ks(p)
{

	inverse=1.0/p->DXM;

	maxF=0.0;
	maxG=0.0; 
	maxK=0.0;
	maxE=0.0;

	sigT=0.9;

	gi=p->W20;
	gj=p->W21;
	gk=p->W22;


}


















