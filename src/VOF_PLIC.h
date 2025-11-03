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
    void vof_transport_COSMIC2D_RK(fdm*,lexer*,ghostcell*,int,int,field&,field&,field&);
    void advectPlane_forCOSMIC2D_RK(fdm*,lexer*,int,int,field&,field&,field&);
    void advectWater_forCOSMIC2D_RK(fdm*,lexer*,int,int,field&,field&,field&);
    void advectPlane_forCOSMIC3D_RK(fdm*,lexer*,int,field&,field&,field&);
    void advectWater_forCOSMIC3D_RK(fdm*,lexer*,int,field&,field&,field&);
    void redistancePhiByPlane_Bonn(fdm*, lexer*);
    double ShortestDistanceOnBoundaryCandidate(fdm*, lexer*, int, int, int, double);
    double ProjectionPointCandidate(fdm*, lexer*, int, int, int, double);
    double IntersectionPointCandidate(fdm*, lexer*, int, int, int, double);
    void symmetric_scheme2D_FCRK3(fdm*, lexer*,ghostcell*, field&, field&, field&);
    void symmetric_scheme3D_FCRK3(fdm*, lexer*,ghostcell*, field&, field&, field&);
    double calcL2vofError2D(fdm*, lexer*, field&, double, double, double, double);
    double calcAlphaFromInput(fdm*, lexer*, double, double, double, double, double, double, double);
    void calcNormalELVIRA2D(fdm*, lexer*, field&);
    void calcNormalMYC2D(fdm*,lexer*, field&);
    void calcNormalMYC3D(fdm*,lexer*, field&);
    void calcNormalMYC3D_V2(fdm*,lexer*, field&);
    void calcNormalMYC2D_V2(fdm*,lexer*, field&);
    void calcNormalMYC2D_V3(fdm*,lexer*, field&);
    int searchMarkerInVicinity(lexer*,fdm*,int,double,int,int,int);
    int searchMarkerAlongDims(lexer*,fdm*,int,double,int,int,int);
    double twoStepVel(lexer*,fdm*,double,double,double);
    
    //COSMICC 3D Subfunctions
    void get_Fx_and_Flux(fdm*,lexer*,ghostcell*,field&);
    void get_Fy_and_Flux(fdm*,lexer*,ghostcell*,field&);
    void get_Fz_and_Flux(fdm*,lexer*,ghostcell*,field&);
    
    void get_Fxy_and_Flux(fdm*,lexer*,ghostcell*,field&);
    void get_Fxz_and_Flux(fdm*,lexer*,ghostcell*,field&);
    void get_Fyx_and_Flux(fdm*,lexer*,ghostcell*,field&);
    void get_Fyz_and_Flux(fdm*,lexer*,ghostcell*,field&);    
    void get_Fzx_and_Flux(fdm*,lexer*,ghostcell*,field&);
    void get_Fzy_and_Flux(fdm*,lexer*,ghostcell*,field&);

    void get_Flux_xyz(fdm*,lexer*,ghostcell*,field&);
    void get_Flux_xzy(fdm*,lexer*,ghostcell*,field&);
    void get_Flux_yxz(fdm*,lexer*,ghostcell*,field&);
    void get_Flux_yzx(fdm*,lexer*,ghostcell*,field&);
    void get_Flux_zxy(fdm*,lexer*,ghostcell*,field&);
    void get_Flux_zyx(fdm*,lexer*,ghostcell*,field&);
    
    void fieldloop_xy(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_xz(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_yx(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_yz(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_zx(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_zy(fdm*,lexer*,ghostcell*,field&,field&,field&);
    
    void fieldloop_xyz(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_xzy(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_yxz(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_yzx(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_zxy(fdm*,lexer*,ghostcell*,field&,field&,field&);
    void fieldloop_zyx(fdm*,lexer*,ghostcell*,field&,field&,field&);


    field4 Vn_p,Vn_m,Vx_p,Vx_m,Vz_p,Vz_m; // 2D Flux Variables
    field4 V_p,V_m;                       // 3D Flux Variables
    
    field4 F_n,F_x,F_y,F_z;         //lvl 1 fields
    int swtch_x, swtch_y, swtch_z;   //lvl 1 switches
    field4 Flux_x,Flux_y,Flux_z;    //lvl 1 Fluxes
    
    field4 F_xy,F_xz,F_yx,F_yz,F_zx,F_zy;                           //lvl 2 fields
    int swtch_xy, swtch_xz, swtch_yx, swtch_yz, swtch_zx, swtch_zy; //lvl 2 switches
    field4 Flux_xy,Flux_xz,Flux_yx,Flux_yz,Flux_zx,Flux_zy;         //lvl 2 Fluxes
    
    field4 Flux_xyz,Flux_xzy,Flux_yxz,Flux_yzx,Flux_zxy,Flux_zyx;       //lvl 3 Fluxes
    int swtch_xyz,swtch_xzy,swtch_yxz,swtch_yzx,swtch_zxy,swtch_zyx;    //lvl 3 switches
    
    field4 phiaux;
    field4 curv;
    fluid_update *pupdate;
    reini *preini;
    interpolation *ipol;
    field4 compressvol;
    field4 phistep;
    
    field4 VoF,vof_rk1,vof_rk2;
    int gcval_vof, gcval_ro, gcval_visc;
    int gcval_phi;

	int gcval_frac;
	double starttime; 
    double ****nxCoeff, ****nyCoeff, ****nzCoeff;
   
	//- Sweep tracker for alternating starting point
	int sSweep;
    int Sweepdim;
    int sweep;
   
    //- Interface normal vector
    field4 nx;
    field4 ny;
    field4 nz;

    double w_thres, a_thres, corr_thres;
    
    //- Plane distance coefficient
    field4 alpha;
    
    int S_S[6][3];
    int S_2D[2][2];
    
    picard_f *ppicard;
    
};
#endif
