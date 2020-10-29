/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2020 Tobias Martin

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

#ifndef mooring_BARDYN_H_
#define mooring_BARDYN_H_

class mooring_barDyn : public mooring
{
public:

	mooring_barDyn(int);
	virtual ~mooring_barDyn();
	
	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void initialize(lexer*, fdm*, ghostcell*);
	virtual void mooringForces(double&, double&, double&);
	
private:

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorXd;
    typedef Eigen::Matrix<double, 3, 3> Matrix3d;
    typedef Eigen::Matrix<double, 1, 3> Vector3d;
    
    typedef vector<vector<double> > MatrixVd;
    typedef vector<vector<int> > MatrixVi;
	

    // Preprocessing
	void ini_parallel(lexer*, fdm*, ghostcell*);
	
	// Runtime
	void startLoop(lexer*, fdm*, ghostcell*, int&);
    void updateField(lexer*, fdm*, ghostcell*, int);
    
    Eigen::VectorXd timeWeight(lexer*);
    
    void updateAcc(lexer*, fdm*, ghostcell*);  
    void updateTopAcc(lexer*);  
    void fillLinSystem(lexer*, fdm*, ghostcell*);
    void fillLinRhs(lexer*, fdm*, ghostcell*);  
    void fillNonLinSystem(lexer*, fdm*, ghostcell*);
    void fillNonLinRhs(lexer*, fdm*, ghostcell*);
    void limitTension();    
    void gravityForce(lexer*);
    void dragForce(lexer*);
    void bottomForce(lexer*);
    void inertiaForce(lexer*);
    
    void print(lexer*);



	// ----------------
	
	// Parallelisation
	int line;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
	// Material constants
	double w, kappa, C1_, C2_, rho_c, d_c;
	
	// Mesh
	double *l0;
    double L; 
    int nl, niK, nbK, nK, nf;
    MatrixXd x0_, x_, xn_, xnn_, xdot_, xdotn_, xdotnn_, xdotdot_, top_xdot_, top_xdotdot_, A_, forces_knot;
    double dt_, dtn_, dtnn_, t_net, t_net_n;
    VectorXd mass_knot, weight_knot, added_mass, B_, T_, T_old, coeffs_;
    double **coupledField, **coupledFieldn, **fb, **K; 
    int *Pi, *Ni;
	
	// Forces
    double Fx,Fy,Fz;
	int **nfK, *nfbK;
	
	// Print
	char name[100];
	ofstream eTout;
	double printtime;
};

#endif
