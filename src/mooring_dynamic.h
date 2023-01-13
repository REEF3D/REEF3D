/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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

#include"mooring.h"
#include"beam.h"
#include"mooring_Catenary.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field5.h"
#include"fieldint5.h"
#include"vec.h"
#include<fstream>
#include<iostream>
#include<vector>
#include <Eigen/Dense>

using namespace std;

#ifndef MOORING_DYN_H_
#define MOORING_DYN_H

class mooring_dynamic : public mooring, public beam
{
public:
	
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    mooring_dynamic(int);
	virtual ~mooring_dynamic();
	
	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void initialize(lexer*, fdm*, ghostcell*);
	virtual void mooringForces(double&, double&, double&);
    
    void setConstantLoads(Matrix3Xd&, Matrix4Xd&, const Matrix3Xd&, const Matrix3Xd&, const Matrix4Xd&, const Matrix4Xd&);
    void setFieldBC(Matrix3Xd&, Matrix3Xd&, Matrix4Xd&, Matrix4Xd&, Matrix4Xd&, Matrix3Xd&, Matrix4Xd&, Matrix3Xd&, double, int);
    
private:

	// Initialisation
	void ini_parallel(lexer*, fdm*, ghostcell*);

	// Runtime
    void updateFields(lexer*, fdm*, ghostcell*);
    void updateFluidVel(lexer*, fdm*, ghostcell*, int);
    void saveMooringPoint(lexer*);

	// ------ 
	
	// Parallelisation
	int line;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
	// Material
	double gamma, A, E, G, L, rho_c, d_c;
	
	// Mesh
    int Ne;

    // Fields
    Matrix3Xd c_moor, c_moor_n, cdot_moor, cdot_moor_n, cdotdot_moor;
    vector<vector<double> > fluid_vel, fluid_vel_n, fluid_acc;

    // Forces
	double Xme_, Yme_, Zme_;
	
    // Time control
	double phi_mooring, t_mooring_n, t_mooring;

    // Boundary conditions
    Eigen::Vector3d fixPoint, a_O;

	// Print
	ofstream eTout;
    
    // Breaking
    bool broken;
    double breakTension, breakTime;
};

#endif
