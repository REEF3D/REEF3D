/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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

#ifndef MOORING_CATENARY_H_
#define MOORING_CATENARY_H_

#include"mooring.h"
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

class mooring_Catenary : public mooring
{
public:
	mooring_Catenary(int);
	virtual ~mooring_Catenary();
	
	virtual void start(lexer*, ghostcell*);
	virtual void initialize(lexer*, ghostcell*);
	virtual void mooringForces(double&, double&, double&);
	
	void getForce(lexer*, ghostcell*, double&, double&);
	void getShape(lexer*, ghostcell*, double*&, double*&, double*&, double*&);
	void iniShape(lexer*, ghostcell*, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&);
	
private:	

    void calcForce(lexer*, ghostcell*);

	// Print
	void print(lexer*);
	void buildLine(lexer*);
	
	// ------ 
	
	// Line number
	int line;
	
	// Material constants
	double w, L, lms, rho_c, EA;
	
	// Mesh
	int H;
	double *x, *y, *z, *B, *F, **A;
	double xs, ys, zs, dx, dy, dz, dxy_aim;
	
	// Forces
	double *T;	
	double FH, FV, FH_0, FV_0, f_Fh, ddf_Fh;
	double Xme_, Yme_, Zme_;
	
	// Print
	char name[100];
	ofstream eTout;
	double printtime;
	
    // Breaking
    bool broken;
    double breakTension, breakTime, curr_time;
};

#endif
