/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fdm_ptf.h"
#include"lexer.h"

fdm_ptf::fdm_ptf(lexer *p):fdm(p),Fifsf(p),K(p),etaloc(p)
{
    
	maxF=0.0;
	maxG=0.0; 
	maxH=0.0;
    
    p->Darray(U_,p->imax*p->jmax*(p->kmax+2));
    p->Darray(V_,p->imax*p->jmax*(p->kmax+2));
    p->Darray(W_,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fi_,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Uin_,p->imax*p->jmax*(p->kmax+2));
    
	gi=p->W20;
	gj=p->W21;
	gk=p->W22;
    
    C4.allocate(p);
    C4a.allocate(p);
    C6.allocate(p);
}

fdm_ptf::~fdm_ptf()
{
    
}