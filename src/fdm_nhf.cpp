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

#include"fdm_nhf.h"
#include"lexer.h"
#include"nhflow_memory.h"

fdm_nhf::fdm_nhf(lexer *p) :  eta(p),etaloc(p),
                              wet_n(p),breaking(p),breaklog(p),bc(p),
                              nodeval2D(p),eta_n(p),WL(p),
                              detadt(p),detadt_n(p),un(p),vn(p),dudt(p),
                              bed(p),depth(p),K(p),
                              Ex(p),Ey(p),Exx(p),Eyy(p),
                              Bx(p),By(p),Bxx(p),Byy(p),
                              hx(p),hy(p),
                              ks(p),
                              coastline(p),PORS(p),vb(p),vbk(p),
                              test2D(p),fs(p),
                              breaking_print(p),Hs(p),
                              rhsvec(p),rvec(p,0),xvec(p,0),N(p,0),M(p),
                              ETAs(p),ETAn(p),ETAe(p),ETAw(p),
                              Ds(p),Dn(p),De(p),Dw(p),dfx(p),dfy(p)
{    
    // largest forcing terms (reset by nhflow_timestep; AMR patch grids start without one)
    maxF=maxG=maxH=0.0;

    p->Darray_lazy(p->sig,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(p->sigx,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(p->sigy,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigz,p->imax*p->jmax);
    p->Darray_lazy(p->sigt,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(p->sigxx,p->imax*p->jmax*(p->kmax+2));
    
    
    p->Darray_lazy(U,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(V,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(W,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(omegaF,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(UH,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(VH,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(WH,p->imax*p->jmax*(p->kmax+2));

    p->Darray_lazy(P,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(RO,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(VISC,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(EV,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(EV0,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(EVH,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(F,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(G,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(H,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(L,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(Fext,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Gext,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Hext,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(POR,p->imax*p->jmax*(p->kmax+2));
    // the fields of some options only (nhflow_memory.h), nullptr otherwise
    if(nhflow_vrans_fields(p))
    p->Darray_lazy(PORPART,p->imax*p->jmax*(p->kmax+2));
    
    if(nhflow_dem_field(p))
    {
    p->Darray_lazy(PORDEM,p->imax*p->jmax*(p->kmax+2));
    for(int n=0; n<p->imax*p->jmax*(p->kmax+2); ++n)
    PORDEM[n] = 1.0;
    }
    
    if(nhflow_test_field(p))
    p->Darray_lazy(test,p->imax*p->jmax*(p->kmax+2));
    
    if(nhflow_kin_field(p))
    p->Darray_lazy(KIN,p->imax*p->jmax*(p->kmax+2));
    
    if(nhflow_conc_field(p))
    p->Darray_lazy(CONC,p->imax*p->jmax*(p->kmax+2));

    // tracers of the semi-implicit NHFLOW (A 600): only when switched on
    if(p->A600>=1)
    p->Darray_lazy(SALT,p->imax*p->jmax*(p->kmax+2));
    if(p->A600>=2)
    p->Darray_lazy(TEMP,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(SOLID,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(FB,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(FHB,p->imax*p->jmax*(p->kmax+2));
    if(nhflow_vrans_fields(p))
    p->Darray_lazy(PORSTRUC,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(Fx,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Fy,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Fz,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(FEx,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(FEy,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(Ss,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Sn,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Se,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Sw,p->imax*p->jmax*(p->kmax+2));
    
    // face states (fdm_nhf.h): the normal velocities, the vertical faces of the velocities and the horizontal
    // faces of the conserved variables; the tangential velocities of the horizontal faces only where they are read
    p->Darray_lazy(Us,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Un,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Ub,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Ut,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(Ve,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Vw,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Vb,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Vt,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(Wb,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Wt,p->imax*p->jmax*(p->kmax+2));
    
    if(nhflow_tangential_faces(p))
    {
    p->Darray_lazy(Ue,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Uw,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Vs,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Vn,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Ws,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Wn,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(We,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(Ww,p->imax*p->jmax*(p->kmax+2));
    }

    p->Darray_lazy(UHs,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(UHn,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(UHe,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(UHw,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(VHs,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(VHn,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(VHe,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(VHw,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray_lazy(WHs,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(WHn,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(WHe,p->imax*p->jmax*(p->kmax+2));
    p->Darray_lazy(WHw,p->imax*p->jmax*(p->kmax+2));
    
    p->Iarray(NODEVAL,p->imax*p->jmax*(p->kmax+3));

}













