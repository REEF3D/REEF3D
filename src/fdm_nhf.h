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

#ifndef FDM_NHF_H_
#define FDM_NHF_H_

#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"
#include"increment.h"
#include"vec.h"
#include"vec2D.h"
#include"matrix_diag.h"
#include"matrix2D.h"

class lexer;
class nhflow_thinbody;
class seastate_nhflow;

using namespace std;

class fdm_nhf : public increment
{
public:

    fdm_nhf(lexer*);
   
    int *NODEVAL;
    
    slice4 eta,eta_n,WL,detadt,detadt_n,un,vn,dudt;
    slice4 bed,depth;
    slice4 K;
    sliceint4 etaloc,wet_n,breaking,breaklog,bc,nodeval2D;
    
    slice4 Ex,Ey;
    slice4 Exx,Eyy;
    slice4 Bx,By;
    slice4 Bxx,Byy;
    
    slice4 hx,hy;
    slice4 ks;
    slice4 coastline;
    slice4 vb;
    slice4 test2D;
    slice4 fs;
    
    slice4 breaking_print,Hs;
    
    // NHFLOW
    
    vec rhsvec;
    vec2D xvec,rvec;
    
    // 3D array
    double *U,*V,*W,*omegaF;
    double *UH,*VH,*WH;
    
    double *P,*RO,*VISC,*EV,*EV0;
    double *Pbc = nullptr;   // iowave: non-hydrostatic pressure of the incoming waves at an open x- edge (Poisson)
    double *F,*G,*H,*L;
    double *Fext,*Gext,*Hext;
    double *POR,*PORPART;
    double *MBETA = nullptr;    // membrane mobility in the pressure Poisson equation (X 330), else unallocated
    double *MCHI = nullptr;     // fraction of the prescribed static pressure below membrane floors (X 330)
    double *MRCX = nullptr;     // membrane (X 330): face correction velocities of the last projection,
    double *MRCY = nullptr;     //                   Rhie-Chow continuity flux next to the membrane
    int MPROJ = 1;              // membrane (X 330): projections per stage (membrane.dat: projections)
    double *MBX = nullptr;      // membrane (X 330), 'mobility link': mobility of the link cell (i,j,k) -> (i+1,j,k),
    double *MBY = nullptr;      //   (i,j,k) -> (i,j+1,k) and of the vertical link node k -> k+1 through cell k;
    double *MBZ = nullptr;      //   unallocated in the default layer mode (mobilities from MBETA)
    nhflow_thinbody *thinbody = nullptr;   // sharp thin bodies (membrane.dat 'mobility sharp'): wall fluxes, projection
                                           // right-hand side and hydrostatic head below closed floors, nhflow_thinbody.h
    int solid_flux = 0;         // 1: no continuity flux through the faces of cells with p->DF < 0 (FEM structures, Z 30)
    double *PORDEM;         // porosity of the REEF3D::DEM particles (E 28), 1 without
    double *test;
    double *KIN;
    double *CONC;
    
    double *SOLID,*FB,*FHB;
    double *PORSTRUC;
    
    double *Fx,*Fy,*Fz;
    double *FEx,*FEy,*FSW,*DWDT;
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
    
    // REEF3D::SEASTATE coupling (A 750 1), nullptr without
    seastate_nhflow *wave = nullptr;
};

#endif
