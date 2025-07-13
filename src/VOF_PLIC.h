/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Tobias Martin, Fabian Knoblauch
--------------------------------------------------------------------*/

#ifndef VOF_PLIC_H_
#define VOF_PLIC_H_

#include"freesurface.h"
#include"gradient.h"
#include"norm_vec.h"
#include"reini.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"interpolation.h"

class picard_f;
class heat;
class fluid_update;

using namespace std;

class VOF_PLIC : public freesurface, gradient, norm_vec
{
public:
	VOF_PLIC(lexer*, fdm*, ghostcell*,heat*);
	virtual ~VOF_PLIC();
	virtual void update(lexer*,fdm*,ghostcell*,field&);
    
    virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
    void RKcalcL(fdm*,lexer*,ghostcell*, field&, field&, field&);
    void RK_redistance(fdm*,lexer*,ghostcell*);
    void updatePhasemarkers(lexer*,fdm*,ghostcell*,field&);
    void updatePhasemarkersCompression(lexer*,fdm*,ghostcell*,field&);
    void updatePhasemarkersCorrection(lexer*,fdm*,ghostcell*,field&);
    void calculateSubFractions(lexer*,fdm*,ghostcell*,field&);
    void surface_tension2D(lexer*,fdm*,ghostcell*,int);
    void updatePlaneData(lexer*,fdm*,ghostcell*,field&);
    double return_alpha_reconstructPlane_alt(fdm*, lexer*,field&,int,int,int);
    
	
private:	
    void iniphi(fdm*, lexer*,ghostcell*);
	void iniphi_io(fdm*, lexer*,ghostcell*);
    void iniphi_box(lexer*,fdm*,ghostcell*);
    void iniphi_surfarea(lexer*,fdm*,ghostcell*);
    int conv(double);
	
	void reconstructPlane(fdm*, lexer*);
	double calcAlpha(fdm*, double&,  double&,  double&);

	void ininorVecLS(lexer*);
	void calcNormalFO(fdm*, lexer*, field&);
	void calcNormalLS(fdm*, lexer*, field&);
    void calcNormalWENO(fdm*, lexer*, field&);
    void calcNormalPhi(fdm*, lexer*);
		
    //Alternative version by Fabian
    void calcNormalWeymouth(fdm*, lexer*, field&);
    void calcNormalWang(fdm*, lexer*);
    void calcNormalMassCentre(fdm*, lexer*, field&);
    void reconstructPlane_alt(fdm*, lexer*,field&);
    double calculateVolume(double,double,double,double,double,double,double);
    void vof_transport_COSMIC2D_RK(fdm*,lexer*,int,int,field&,field&,field&);
    void simpleNormal_Bonn(fdm*, lexer*);
    void advectPlane_forCOSMIC2D_RK(fdm*,lexer*,int,int,field&,field&,field&);
    void advectWater_forCOSMIC2D_RK(fdm*,lexer*,int,int,field&,field&,field&);
    void redistancePhiByPlane_Bonn(fdm*, lexer*);
    double ShortestDistanceOnBoundaryCandidate(fdm*, lexer*, int, int, int, double);
    double ProjectionPointCandidate(fdm*, lexer*, int, int, int, double);
    double IntersectionPointCandidate(fdm*, lexer*, int, int, int, double);
    void symmetric_scheme2D_FCRK3(fdm*, lexer*,ghostcell*, field&, field&, field&);
    double calcL2vofError2D(fdm*, lexer*, field&, double, double, double, double);
    double calcAlphaFromInput(fdm*, lexer*, double, double, double, double, double, double, double);
    void calcNormalELVIRA2D(fdm*, lexer*, field&);
    void calcNormalMYC2D(fdm*,lexer*, field&);
    void calcNormalMYC2D_V2(fdm*,lexer*, field&);
    void calcNormalMYC2D_V3(fdm*,lexer*, field&);
    int searchMarkerInVicinity(lexer*,fdm*,int,double,int,int,int);
    int searchMarkerAlongDims(lexer*,fdm*,int,double,int,int,int);
    double twoStepVel(lexer*,fdm*,double,double,double);
    
   
    field4 V_w_p;
    field4 V_w_m;
    field4 Vx_p; // Vx_+1/2
    field4 Vx_m; // Vx_-1/2
    field4 Vz_p; // Vz_+1/2
    field4 Vz_m; // Vz_-1/2
    field4 Vn_p; // Vn_+1/2
    field4 Vn_m; // Vn_-1/2
    field4 F_x; // Fx in vof transport scheme
    field4 F_z; // Fz in vof transport scheme
    field4 F_n; // old vof field
    field4 F_new; // F_n+1 new vof field
    field4 Flux_x; // delta_t*d/dx(uF_n)
    field4 Flux_z; // delta_t*d/dz(wF_n)
    field4 Crossflux_xz; // delta_t*d/dz(wFx)
    field4 Crossflux_zx; // delta_t*d/dx(uFz)
    field4 phival;
    field4 alphastore;
    field4 phistep;
    field4 phiS0;
    field4 phiS1;
    field4 phiS2;
    field4 vofstep;
    field4 vofS0;
    field4 vofS1;
    field4 vofS2;
    field4 phiaux;
    field4 curv;
    fluid_update *pupdate;
    reini *preini;
    interpolation *ipol;
    field4 compressvol;
    
    field4 VoF,vof_rk1,vof_rk2;
    int gcval_vof, gcval_ro, gcval_visc;
    int gcval_phi;

	int gcval_frac;
	double starttime; 
   
	//- Sweep tracker for alternating starting point
	int sSweep;
    int Sweepdim;
    int sweep;
   
    //- Interface normal vector
    field4 nx;
    field4 ny;
    field4 nz;
    
	double ****nxCoeff, ****nyCoeff, ****nzCoeff;
    
    double w_thres, a_thres, corr_thres;
    
    //- Plane distance coefficient
    field4 alpha;
    
    //- Volume fractions in the cell and its neighbours
    field4 vof1;
    field4 vof2;
    field4 vof3;
    
    int S_S[6][3];
    int S_2D[2][2];
    
    picard_f *ppicard;
    
};
#endif
