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

#include"fdm_nhf.h"
#include"lexer.h"

fdm_nhf::fdm_nhf(lexer *p) :  eta(p),etaloc(p),
                              wet_n(p),breaking(p),breaklog(p),bc(p),
                              nodeval2D(p),
                              eta_n(p),WL(p),detadt(p),detadt_n(p),
                              bed(p),depth(p),K(p),
                              Ex(p),Ey(p),Exx(p),Eyy(p),
                              Bx(p),By(p),Bxx(p),Byy(p),
                              hx(p),hy(p),
                              coastline(p),vb(p),
                              test2D(p),fs(p),
                              breaking_print(p),Hs(p),
                              rhsvec(p),rvec(p),xvec(p),N(p),M(p),
                              ETAs(p),ETAn(p),ETAe(p),ETAw(p),
                              Ds(p),Dn(p),De(p),Dw(p),dfx(p),dfy(p)
{    
    p->Darray(p->sig, p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigz,p->imax*p->jmax);
    p->Darray(p->sigt,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+2));
    
    
    p->Darray(U,p->imax*p->jmax*(p->kmax+2));
    p->Darray(V,p->imax*p->jmax*(p->kmax+2));
    p->Darray(W,p->imax*p->jmax*(p->kmax+2));
    p->Darray(omegaF,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(UH,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VH,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WH,p->imax*p->jmax*(p->kmax+2));

    p->Darray(P,p->imax*p->jmax*(p->kmax+2));
    p->Darray(RO,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VISC,p->imax*p->jmax*(p->kmax+2));
    p->Darray(EV,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(F,p->imax*p->jmax*(p->kmax+2));
    p->Darray(G,p->imax*p->jmax*(p->kmax+2));
    p->Darray(H,p->imax*p->jmax*(p->kmax+2));
    p->Darray(L,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(porosity,p->imax*p->jmax*(p->kmax+2));
    p->Darray(test,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Fx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fy,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fz,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Fs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fw,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Ss,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Sn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Se,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Sw,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(SSx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(SSy,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(Us,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Un,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ue,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Uw,p->imax*p->jmax*(p->kmax+2));  
    p->Darray(Ub,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ut,p->imax*p->jmax*(p->kmax+2));  
    
    p->Darray(Vs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ve,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vw,p->imax*p->jmax*(p->kmax+2));  
    p->Darray(Vb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Vt,p->imax*p->jmax*(p->kmax+2));  
    
    p->Darray(Ws,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Wn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(We,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Ww,p->imax*p->jmax*(p->kmax+2));  
    p->Darray(Wb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Wt,p->imax*p->jmax*(p->kmax+2)); 


    p->Darray(UHs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(UHn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(UHe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(UHw,p->imax*p->jmax*(p->kmax+2));  
    p->Darray(UHb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(UHt,p->imax*p->jmax*(p->kmax+2));  
    
    p->Darray(VHs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHw,p->imax*p->jmax*(p->kmax+2));  
    p->Darray(VHb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VHt,p->imax*p->jmax*(p->kmax+2));  
    
    p->Darray(WHs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHw,p->imax*p->jmax*(p->kmax+2));  
    p->Darray(WHb,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WHt,p->imax*p->jmax*(p->kmax+2)); 
    
    p->Iarray(NODEVAL,p->imax*p->jmax*(p->kmax+3));

}













