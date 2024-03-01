/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"
#include"sliceint5.h"
#include"increment.h"
#include"vec.h"
#include"vec2D.h"
#include"matrix_diag.h"
#include"matrix2D.h"
#include"cpt2D.h"

class lexer;

#ifndef FDM_NHF_H_
#define FDM_NHF_H_

using namespace std;

class fdm_nhf : public increment
{
public:

    fdm_nhf(lexer*);
   
    int *NODEVAL;
    
    slice4 eta,eta_n,WL,detadt,detadt_n;
    slice4 bed,depth;
    slice4 K;
    sliceint4 etaloc,wet_n,breaking,breaklog,bc,nodeval2D;
    
    slice4 Ex,Ey;
    slice4 Exx,Eyy;
    slice4 Bx,By;
    slice4 Bxx,Byy;
    
    slice4 hx,hy;
    slice4 coastline;
    slice4 vb;
    slice4 test2D;
    
    slice4 breaking_print,Hs;
    
    // NHFLOW
    
    vec rhsvec;
    vec2D xvec,rvec;
    
    // 3D array
    double *U,*V,*W,*omegaF;
    double *UH,*VH,*WH;
    
    double *P,*RO,*VISC,*EV;
    double *F,*G,*H,*L;
    double *porosity;
    double *test;
    
    double *Fx,*Fy,*Fz;
    double *Fs,*Fn,*Fe,*Fw;
    double *Ss,*Sn,*Se,*Sw;
    double *SSx,*SSy;
    
    double *Un,*Us,*Ue,*Uw,*Ub,*Ut;
    double *Vn,*Vs,*Ve,*Vw,*Vb,*Vt;
    double *Wn,*Ws,*We,*Ww,*Wb,*Wt;
    
    double *UHn,*UHs,*UHe,*UHw,*UHb,*UHt;
    double *VHn,*VHs,*VHe,*VHw,*VHb,*VHt;
    double *WHn,*WHs,*WHe,*WHw,*WHb,*WHt;
    
    slice1 ETAs,ETAn;
    slice2 ETAe,ETAw;
    slice1 Ds,Dn;
    slice2 De,Dw;
    slice1 dfx;
    slice2 dfy;

    matrix2D N;
	matrix_diag M;    
    
    double gi,gj,gk;
    double maxF,maxG,maxH;
    double wd_criterion;
    
};

#endif
