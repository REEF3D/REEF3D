/*--------------------------------------------------------------------
REEF3D
Copyright 2019 Tobias Martin

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

#include"net.h"

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

#ifndef NET_QUASISTATIC_H_
#define NET_QUASISTATIC_H_

class net_QuasiStatic : public net
{
public:
	net_QuasiStatic(int);
	virtual ~net_QuasiStatic();
    
	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void initialize(lexer*, fdm*, ghostcell*);
	virtual void netForces(double&, double&, double&);
    
    
private:

    // Preprocessing
    void genericNet();
    void iniInnerKnots();
    void iniBoundaryKnots();  
    void stretch();
    void iniLSE(lexer*);
	void ini_parallel(lexer*, fdm*, ghostcell*);
    
	// Runtime
    void solveGauss(lexer*, double**&, double**&, double**&);
    void updateVel(lexer*, fdm*, ghostcell*, int);
	void getC(double, double*&);
	void updateLength();
    void print(lexer*);
    void buildNet(lexer*);
	
	// ------ 
	
	// Parallelisation
	int nNet;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
	// Material constants
	double EA, w, rho_c, d_c;
	
	// Mesh
	double origin_x, origin_y, origin_z, phi, theta, psi;
	double *l0, *l;
    double L, b, lm, Fg, beta, gamma; 
    int n, m, niK, nK, nf;
    double **fi, **fb, **K, **A, **B, **Bh, **A_, **B_, **K_;
    int *Pb, *Nb, *Pi, *Ni;
	
	// Forces
	double **v, **c, **e, **e_q, **e_d, **e_l;
	int **nfK;

	// Print
	double *x, *y, *z, *T;
	char name[100];
	ofstream eTout;
	double printtime;
};

#endif
