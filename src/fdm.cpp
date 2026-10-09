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
            ro(p),visc(p),
            phi(p),
            vof(p),vof_nt(p,p->F80>0),vof_nb(p,p->F80>0),vof_st(p,p->F80>0),vof_sb(p,p->F80>0),phasemarker(p,p->F80>0),
            vof_nte(p,p->F80>0),vof_ntw(p,p->F80>0),vof_nbe(p,p->F80>0),vof_nbw(p,p->F80>0),vof_ste(p,p->F80>0),vof_stw(p,p->F80>0),vof_sbe(p,p->F80>0),vof_sbw(p,p->F80>0),
            conc(p),
            topo(p),solid(p),
            test(p),
            fb(p),fbh1(p),fbh2(p),fbh3(p),fbh4(p),fbh5(p),
            porosity(p),porpart(p),porA(p),porB(p),
            walld(p),
            nodeval(p),nodeval2D(p),etaloc(p),
            eta(p),eta_n(p),depth(p),WL(p),
            Fifsf(p),K(p),
            P(p),Q(p),bed(p),
            rhsvec(p),M(p),
            nX(p,p->F80>0),nY(p,p->F80>0),nZ(p,p->F80>0),Alpha(p,p->F80>0)
            
{
	maxF=0.0;
	maxG=0.0; 
	maxH=0.0;
    
	gi=p->W20;
	gj=p->W21;
	gk=p->W22;
}


















