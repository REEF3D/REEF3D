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

#ifndef mooring_BARQUASISTATIC_H_
#define mooring_BARQUASISTATIC_H_

#include"mooring.h"
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


using namespace std;

class mooring_barQuasiStatic : public mooring
{
public:
	mooring_barQuasiStatic(int);
	virtual ~mooring_barQuasiStatic();
	
	virtual void start(lexer*, ghostcell*);
	virtual void initialize(lexer*, ghostcell*);
	virtual void mooringForces(double&, double&, double&);
	
private:

	// Preprocessing
	void ini_parallel(lexer*, ghostcell*);
	
	// Runtime
	vector< vector<double> > solveGauss(vector< vector<double> >, vector< vector<double> >);
    void updateVel(lexer*, ghostcell*, int);
 	void getVel(lexer*, ghostcell*);
        vector<double> getC(double);
    void print(lexer*,ghostcell*);
    void buildLine(lexer*,ghostcell*);
    void checkBottom(lexer*,ghostcell*);
	
	// ------ 
	
	// Parallelisation
	int line;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
	// Material constants
	double w, EA, L, rho_c, d_c;
	
	// Mesh
	int sigma;
	double dx, dy, dz;
	vector<double> l0, l;
	double *x, *y, *z;
	
	// Fields
        vector< vector<double> > u_, v_, w_, A, B, f, R, v;
	
	// Forces
	double *T, *Fb;
	vector< vector<double> > e, c_coeff, e_q, e_d, e_l;
	double Xme_, Yme_, Zme_;
	
	// Print
	char name[100];
	ofstream eTout;
	double printtime;

    // Catenary
	mooring_Catenary *pcatenary;

    // Breaking
    bool broken;
    double breakTension, breakTime, curr_time;
};

#endif
