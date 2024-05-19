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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"freesurface.h"
#include"gradient.h"
#include"norm_vec.h"
#include"reini.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

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
	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
	virtual void update(lexer*,fdm*,ghostcell*,field&);
	
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
	
	
    fluid_update *pupdate;
    reini *reini_;

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
};
#endif
