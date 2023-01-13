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

#include"fdm_nhf.h"
#include"lexer.h"

fdm_nhf::fdm_nhf(lexer *p) :  press(p),nodeval(p),eta(p),etaloc(p),
                              wet_n(p),breaking(p),breaklog(p),bc(p),
                              eta_n(p),WL(p),WL_n(p),
                              bed(p),depth(p),K(p),
                              Fx(p),Fy(p),
                              Ex(p),Ey(p),Exx(p),Eyy(p),
                              Bx(p),By(p),Bxx(p),Byy(p),
                              hx(p),hy(p),
                              wbed(p),dwdt(p),
                              coastline(p),vb(p),
                              nodeval2D(p),breaking_print(p),
                              rhsvec(p),rvec(p),xvec(p),N(p),M(p)
{    
    p->Darray(U,p->imax*p->jmax*(p->kmax+2));
    p->Darray(V,p->imax*p->jmax*(p->kmax+2));
    p->Darray(W,p->imax*p->jmax*(p->kmax+2));
    p->Darray(omega,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(P,p->imax*p->jmax*(p->kmax+2));
    p->Darray(ro,p->imax*p->jmax*(p->kmax+2));
    p->Darray(visc,p->imax*p->jmax*(p->kmax+2));
    p->Darray(eddyv,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(F,p->imax*p->jmax*(p->kmax+2));
    p->Darray(G,p->imax*p->jmax*(p->kmax+2));
    p->Darray(H,p->imax*p->jmax*(p->kmax+2));
    p->Darray(L,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(porosity,p->imax*p->jmax*(p->kmax+2));
    p->Darray(test,p->imax*p->jmax*(p->kmax+2));

    C4.allocate(p);
}













