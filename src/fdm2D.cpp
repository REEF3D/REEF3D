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

#include"fdm2D.h"
#include"lexer.h"
#include"sflow_memory.h"

fdm2D::fdm2D(lexer *p)
			: U(p),V(p),W(p),UH(p),VH(p),WH(p),UA(p),VA(p),MX(p),MY(p),
            Fx(p),Fs(p,sflow_face_fluxes(p)),Fn(p,sflow_face_fluxes(p)),FEx(p),
            Fy(p),Fe(p,sflow_face_fluxes(p)),Fw(p,sflow_face_fluxes(p)),FEy(p),
    
            Ss(p),Sn(p),
            Se(p),Sw(p),
    
            Un(p),Us(p),
            Ue(p,sflow_tangential_faces(p)),Uw(p,sflow_tangential_faces(p)),
    
            Vn(p,sflow_tangential_faces(p)),Vs(p,sflow_tangential_faces(p)),
            Ve(p),Vw(p),
    
            Wn(p,sflow_w_faces(p)),Ws(p,sflow_w_faces(p)),
            We(p,sflow_w_faces(p)),Ww(p,sflow_w_faces(p)),
    
            UHn(p),UHs(p),
            UHe(p),UHw(p),
    
            VHn(p),VHs(p),
            VHe(p),VHw(p),
    
            WHn(p,sflow_wh_faces(p)),WHs(p,sflow_wh_faces(p)),
            WHe(p,sflow_wh_faces(p)),WHw(p,sflow_wh_faces(p)), 
            
            QUn(p,sflow_bous_faces(p)),QUs(p,sflow_bous_faces(p)),
            QUe(p,sflow_bous_faces(p)),QUw(p,sflow_bous_faces(p)),
            QVn(p,sflow_bous_faces(p)),QVs(p,sflow_bous_faces(p)),
            QVe(p,sflow_bous_faces(p)),QVw(p,sflow_bous_faces(p)),
            
            WL(p),
            ETAs(p),ETAn(p),
            ETAe(p),ETAw(p),
            
            Ds(p),Dn(p),
            De(p),Dw(p),
            
            dfx(p),dfy(p),
    
    
            eta(p),eta_n(p),
            P(p),Q(p),
            F(p),G(p),H(p),L(p),
            Fext(p),Gext(p),Hext(p),
            ws(p),
            press(p),
            eddyv(p),kin(p),
            bed(p),bed0(p),depth(p),
            solidbed(p),topobed(p),
            bednode(p),
            hx(p),hy(p),hp(p),
            xvec(p),rhsvec(p),M(p),
            test(p),pbh(p),Hs(p),fs(p),
            fx(p),fy(p),fz(p),
            breaking(p),breaking_print(p),vb(p),
            wet1(p),deep1(p),wet2(p),deep2(p),
            nodeval(p),
            cmu(0.09),
            ks(p),spiral(p)
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


















