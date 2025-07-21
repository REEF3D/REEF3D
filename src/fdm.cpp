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

#include"fdm.h"
#include"lexer.h"

fdm::fdm(lexer *p) :
            u(p),F(p),Fext(p),
            v(p),G(p),Gext(p),
            w(p),H(p),Hext(p),
            press(p),
            Fi(p),
            eddyv(p),
            L(p),
            ro(p),dro(p),visc(p),
            phi(p),
            vof(p),vof_nt(p),vof_nb(p),vof_st(p),vof_sb(p),phasemarker(p),
            vof_nte(p),vof_ntw(p),vof_nbe(p),vof_nbw(p),vof_ste(p),vof_stw(p),vof_sbe(p),vof_sbw(p),
            conc(p),
            topo(p),solid(p),
            test(p),
            fb(p),fbh1(p),fbh2(p),fbh3(p),fbh4(p),fbh5(p),
            porosity(p),porpart(p),
            walld(p),
            nodeval(p),nodeval2D(p),etaloc(p),
            eta(p),eta_n(p),depth(p),
            Fifsf(p),K(p),
            P(p),Q(p),bed(p),
            rhsvec(p),M(p),
            nX(p),nY(p),nZ(p),Alpha(p)
            
{
	maxF=0.0;
	maxG=0.0; 
	maxH=0.0;
    
	gi=p->W20;
	gj=p->W21;
	gk=p->W22;
}


















