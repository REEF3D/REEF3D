/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2025 Tobias Martin

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#ifndef NET_BARDYN_H_
#define NET_BARDYN_H_

#include"net.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4a.h"
#include"fieldint4.h"
#include"field5.h"
#include"fieldint5.h"
#include"vec.h"
#include<fstream>
#include<iostream>
#include<vector>
#include"boundarycheck.h"

class reinidisc;

using namespace std;

class net_barDyn : public net, public boundarycheck
{
public:
	net_barDyn(int, lexer*);
	virtual ~net_barDyn();
    
	virtual void start_cfd(lexer*, fdm*, ghostcell*, double, Eigen::Matrix3d);
    virtual void start_nhflow(lexer*, fdm_nhf*, ghostcell*, double, Eigen::Matrix3d);
    
	virtual void initialize_cfd(lexer*, fdm*, ghostcell*);
    virtual void initialize_nhflow(lexer*, fdm_nhf*, ghostcell*);
	virtual void netForces(lexer*, double&, double&, double&, double&, double&, double&);
    
    virtual const EigenMat& getLagrangePoints(){return lagrangePoints;} 
    virtual const EigenMat& getLagrangeForces(){return lagrangeForces;} 
    virtual const EigenMat& getCollarVel(){return collarVel;} 
    virtual const EigenMat& getCollarPoints(){return collarPoints;} 

    
private:
    
    // -------------------------------
	// Runtime
	void startLoop(lexer*, ghostcell*, int&);
    void update_velocity_cfd(lexer*, fdm*, ghostcell*);
    void update_velocity_nhflow(lexer*, fdm_nhf*, ghostcell*);
    
    void updateField_cfd(lexer*, fdm*, ghostcell*, int);
    void updateField_nhflow(lexer*, fdm_nhf*, ghostcell*, int);
    
    void coupling_dlm_cfd(lexer*, fdm*, ghostcell*);
    void coupling_dlm_nhflow(lexer*, fdm_nhf*, ghostcell*);
    
    // -------------------------------
    
    // Preprocessing
    void cone_ini(lexer*, ghostcell*); 
    void cyl_ini(lexer*, ghostcell*); 
    void wall_ini(lexer*, ghostcell*); 
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorXd;
    typedef Eigen::Matrix<double, 3, 3> Matrix3d;
    typedef Eigen::Matrix<double, 1, 3> Vector3d;
    
    typedef vector<vector<double> > MatrixVd;
    typedef vector<vector<int> > MatrixVi;
    
    void print(lexer*);
    
	void buildNet_cyl(lexer*);
    void buildNet_wall(lexer*);
    
    Eigen::VectorXd timeWeight(lexer*);
    
    void updateAcc(lexer*, ghostcell*);  
    void updateTopAcc(lexer*);  
    void fillLinSystem(lexer*, ghostcell*);
    void fillLinRhs(lexer*, ghostcell*);  
    void fillNonLinSystem(lexer*, ghostcell*);
    void fillNonLinRhs(lexer*, ghostcell*);
    void limitTension();    
    
    void getForces(lexer*);
    void gravityForce(lexer*);
    void dragForce(lexer*);
    void inertiaForce(lexer*);
    Eigen::Vector3d screenForce(lexer*, const double&, const Vector3d&, const Vector3d&, const double&, const int, const int);
    void screenForceCoeff(lexer*,double&, double&, const double&, const double&, const double&);
    
    
    void triangulation(lexer*, ghostcell*);
    void create_triangle
        (
            MatrixVd&, MatrixVd&, MatrixVd&,
            const double&,const double&,const double&,
            const double&,const double&,const double&,const double&,
            const double&,const double&,const double&,const double&,
            const double&
        );
    
    
	// ------ 
	
    
	// Parallelisation
	int nNet;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
	// Material constants
	double EA, w, rho_c, l_c, d_c, kappa, C1_, C2_;
	
	// Mesh
	double origin_x, origin_y, origin_z, phi, theta, psi;
	double *l0;
    double L, D, al, ad, beta, gamma; 
    int nd, nl, niK, nbK, nK, nf;
    MatrixXd x0_, x_, xn_, xnn_, xdot_, xdotn_, xdotnn_, xdotdot_, top_xdot_, top_xdotdot_, A_, forces_knot;
    double dt_, dtn_, dtnn_, t_net, t_net_n;
    VectorXd mass_knot, weight_knot, sinker_knot, added_mass, B_, T_, T_old, T_backup, coeffs_;
    double **coupledField, **coupledFieldn, **fb, **K; 
    int *Pb, *Nb, *Pi, *Ni;
    vector<vector<int> > meshID;
    
    // Net mesh
    int tend;  
    MatrixVd tri_x, tri_y, tri_z, tri_vel, tri_forces;
    vector<Eigen::Vector3d> lagrangePoints;    
    vector<Eigen::Vector3d> lagrangeForces;  
    vector<Eigen::Vector3d> collarVel;    
    vector<Eigen::Vector3d> collarPoints;    
	
    // Forces
    double Tne,Fx,Fy,Fz;
	int **nfK, *nfbK;

    // Sinker
    double sinker_m, sinker_l, sinker_d;

    // Knots
    double knot_d;

    // Probe Points
    Eigen::VectorXi probeKnot;

	// Print
	double *x, *y, *z, *T;
	char name[100];
	double printtime;
};

#endif
