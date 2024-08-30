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
Author: Tobias Martin, Fabian Knoblauch
--------------------------------------------------------------------*/

#include"freesurface.h"
#include"gradient.h"
#include"norm_vec.h"
#include"reini.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"interpolation.h"

class picard;
class heat;
class fluid_update;

using namespace std;

#ifndef VOF_PLIC_H_
#define VOF_PLIC_H_

class VOF_PLIC : public freesurface, gradient, norm_vec
{
public:
	VOF_PLIC(lexer*, fdm*, ghostcell*,heat*);
	virtual ~VOF_PLIC();
	virtual void start_old(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
    virtual void start_work(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
	virtual void update(lexer*,fdm*,ghostcell*,field&);
    
    virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
	
private:	
    void iniphi(fdm*, lexer*,ghostcell*);
	void iniphi_io(fdm*, lexer*,ghostcell*);
    void iniphi_box(lexer*,fdm*,ghostcell*);
    void iniphi_surfarea(lexer*,fdm*,ghostcell*);
    int conv(double);
	
	void reconstructPlane(fdm*, lexer*);
	double calcAlpha(fdm*, double&,  double&,  double&);

	void ininorVecLS(lexer*);
	void calcNormalFO(fdm*, lexer*);
	void calcNormalLS(fdm*, lexer*);
    void calcNormalWENO(fdm*, lexer*);
    void calcNormalPhi(fdm*, lexer*);
	
	void advectPlane(fdm*, lexer*, double, double, int);
	void calcFlux(fdm*, lexer*, double&, double&, int);
	
	void updateVolumeFraction(fdm*, lexer*, double, double, int);
	void updateVOF(fdm*, lexer*, int);
	double calcV(const double&, const double&, const double&, const double&, double, double);
	double calcV2(lexer*);
	double F3(const double&);
	
    void redistance(fdm*, lexer*, convection*, ghostcell*,ioflow*,int);
    int calcBoundaryPoint(fdm*, lexer*, int, int, int, field4&);
    int calcProjectionPoint(fdm*, lexer*, double&, double&, double&, int, int, int, field4&);
	void calcSegmentPoint(fdm*, lexer*, double, double, double, int, int, int, field4&);
    double calcDistance(double, double, double, double);	
	
    //Alternative version by Fabian
    void calculateNormal_alt(fdm*, lexer*);
    void reconstructPlane_alt(fdm*, lexer*);
    void advectPlane_alt(fdm*, lexer*,int);
    void advectPlane_altFlux(fdm*, lexer*,double,double,int);
    double calculateVolume(double,double,double,double,double,double,double);
    void updateVOF_alt(fdm*, lexer*,int);
    void advectPlane_sweepless(fdm*, lexer*);
    void updateVOF_sweepless(fdm*, lexer*);
    void advectWater_sweepless(fdm*, lexer*);
    void advectPlane_Weymouth(fdm*, lexer*, int);
    void advectWater_Weymouth(fdm*, lexer*, int);
    void advectPlane_MACHO2D(fdm*, lexer*, int);
    void advectWater_MACHO2D(fdm*, lexer*, int);
    void updateVOF_MACHO2D(fdm*, lexer*, int, int);
    void updateVOF_Weymouth(fdm*, lexer*, int);
    void advectWater_WeymouthNoS(fdm*, lexer*);
    void advectPlane_Wang(fdm*, lexer*, int);
    void transportPhi_Bonn(fdm*,lexer*,int,int);
    void transportVOF_Bonn(fdm*,lexer*,int,int);
    void simpleNormal_Bonn(fdm*, lexer*);
    void advectPlane_forBonnScheme(fdm*, lexer*,int);
    void advectWater_forBonnScheme(fdm*, lexer*,int);
    void redistancePhiByPlane_Bonn(fdm*, lexer*);
    double ShortestDistanceOnBoundaryCandidate(fdm*, lexer*, int, int, int);
    double ProjectionPointCandidate(fdm*, lexer*, int, int, int);
    double IntersectionPointCandidate(fdm*, lexer*, int, int, int);
    
    field4 vof_old;
    field4 V_w_p;
    field4 V_w_m;
    field4 V_w_update;
    field4 V_a_update;
    field4 phival;
    field4 Watersafe;
    field4 V_w_p_star;
    field4 V_w_m_star;
    field4 vof_prevstep;
    field4 V_w_old;
    field4 V_a_old;
    field4 FX_p;
    field4 FX_m;
    field4 FZ_p;
    field4 FZ_m;
    field4 alphastore;
    field4 nxstore;
    field4 nystore;
    field4 nzstore;
    field4 phistep;
    field4 phiS0;
    field4 phiS1;
    field4 phiS2;
    field4 vofstep;
    field4 vofS0;
    field4 vofS1;
    field4 vofS2;
    field4 phiaux;
    fluid_update *pupdate;
    reini *reini_;
    interpolation *ipol;

	int gcval_frac;
	double starttime; 
   
	//- Sweep tracker for alternating starting point
	int sSweep;
   
    //- Interface normal vector
    field4 nx;
    field4 ny;
    field4 nz;
    
	double ****nxCoeff, ****nyCoeff, ****nzCoeff;
    
    //- Plane distance coefficient
    field4 alpha;
    
    //- Volume fractions in the cell and its neighbours
    field4 vof1;
    field4 vof2;
    field4 vof3;
    
    int S_S[6][3];
    
};
#endif
