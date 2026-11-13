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
#include"sliceint5.h"
#include"increment.h"
#include"vec.h"
#include"vec2D.h"
#include"matrix_diag.h"
#include"matrix2D.h"

class lexer;
class nhflow_thinbody;
class nhflow_fsf_hull;
class nhflow_hullbody;
class seastate_nhflow;
class met_forcing;

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
    slice4 ks;
    slice4 coastline;
    slice4 PORS;                // VRANS: porosity at the free surface n(eta), storage of the continuity
    slice4 vb;
    slice4 vbk;     // A 550 4: wavenumber of the breaking crest for the depth profile of vb
    slice4 test2D;
    slice4 fs;
    
    slice4 breaking_print,Hs;
    
    // NHFLOW
    
    vec rhsvec;
    vec2D xvec,rvec;    // not used by NHFLOW: empty (they were 2D vectors of the 3D row count)
    
    // 3D array
    double *U,*V,*W,*omegaF;
    double *UH,*VH,*WH;
    
    double *P,*RO,*VISC,*EV,*EV0;
    double *EVH;   // horizontal eddy viscosity of the HLES split (A 561 > 0), the momentum x/y stresses
    double *Pbc = nullptr;   // iowave: non-hydrostatic pressure of the incoming waves at an open x- edge (Poisson)
    double *F,*G,*H,*L;
    double *Fext,*Gext,*Hext;
    double *POR;
    double *PORPART = nullptr;  // VRANS porous structures (B 200 1)
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
    slice4 *zeta = nullptr;     // X 13 2: level of the hydrostatic pressure terms (the piezometric head under a floating
                                //         body, nhflow_fsf_hull), eta where nullptr
    nhflow_fsf_hull *hull = nullptr;   // X 13 2: the hull as the top of the water columns (nhflow_fsf_hull.h)
    nhflow_hullbody *hullbody = nullptr;   // X 15 8: the hull as a sharp thin body (nhflow_hullbody.h)
    int solid_flux = 0;         // 1: no continuity flux through the faces of cells with p->DF < 0 (FEM structures, Z 30); 2: only p->DF == -1 (direct forcing solids, A 521 0)
    // the fields with nullptr below are allocated only for the options that read them (nhflow_memory.h)
    double *PORDEM = nullptr;   // porosity of the REEF3D::DEM particles (E 10 > 0)
    double *test = nullptr;     // output P 23, work array of A 520 3
    double *KIN = nullptr;      // turbulence models (A 560 > 0)
    double *CONC = nullptr;     // suspended sediment (S 12 > 0), output P 26
    double *SALT = nullptr;     // NHFLOW semi-implicit tracers (A 600 1/2): salinity (psu), cell centres
    double *TEMP = nullptr;     // A 600 2: temperature (deg C)
    
    double *SOLID,*FB,*FHB;
    double *PORSTRUC = nullptr; // VRANS porous structures (B 200 1)
    
    double *Fx,*Fy,*Fz;
    double *FEx,*FEy;
    double *Ss,*Sn,*Se,*Sw;
    
    // face states: s/n the left/right state of the x face i+1/2, e/w of the y face j+1/2, b/t of the z face k+1/2.
    // The tangential velocities of the horizontal faces (Vs..Ww, nullptr) only for HLLC, the curvilinear grid
    // and the semi-implicit momentum (nhflow_tangential_faces); the conserved variables have no z faces.
    double *Un,*Us,*Ub,*Ut;
    double *Ve,*Vw,*Vb,*Vt;
    double *Wb,*Wt;
    double *Vn = nullptr, *Vs = nullptr, *Wn = nullptr, *Ws = nullptr;
    double *Ue = nullptr, *Uw = nullptr, *We = nullptr, *Ww = nullptr;
    
    double *UHn,*UHs,*UHe,*UHw;
    double *VHn,*VHs,*VHe,*VHw;
    double *WHn,*WHs,*WHe,*WHw;
    
    slice1 ETAs,ETAn;
    slice2 ETAe,ETAw;
    slice1 Ds,Dn;
    slice2 De,Dw;
    slice1 dfx;
    slice2 dfy;

    matrix2D N;         // not used by NHFLOW: empty
	matrix_diag M;    
    
    double gi,gj,gk;
    double maxF,maxG,maxH;
    double wd_criterion;
    
    // REEF3D::SEASTATE coupling (A 750 1), nullptr without
    seastate_nhflow *wave = nullptr;

    // REEF3D atmospheric forcing (W 200) of this grid, owned by the momentum class (output)
    met_forcing *met = nullptr;
};

#endif
