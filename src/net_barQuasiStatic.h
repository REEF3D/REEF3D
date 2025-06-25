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

#ifndef NET_BARQUASISTATIC_H_
#define NET_BARQUASISTATIC_H_

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

class net_barQuasiStatic : public net, public boundarycheck
{
public:
	net_barQuasiStatic(int, lexer*);
	virtual ~net_barQuasiStatic();
    
	virtual void start_cfd(lexer*, fdm*, ghostcell*, double,Eigen::Matrix3d)=0;
    virtual void start_nhflow(lexer*, fdm_nhf*, ghostcell*, double,Eigen::Matrix3d)=0;
    
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
    void bag_ini(lexer*, ghostcell*);
    void cyl_ini(lexer*, ghostcell*); 
    void wall_ini(lexer*, ghostcell*); 
    void genericNet();
    void iniInnerKnots();
    void iniBoundaryKnots();  
    void stretch();
    void iniLSE(lexer*);
	void ini_parallel(lexer*, ghostcell*);
    
    
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    typedef Eigen::Matrix<double, 3, 3> Matrix3d;
    typedef Eigen::Matrix<double, 1, 3> Vector3d;
    
    typedef vector<vector<double> > MatrixVd;
    typedef vector<vector<int> > MatrixVi;
    
    
    vector<Eigen::Vector3d> lagrangePoints;    
    vector<Eigen::Vector3d> lagrangeForces;    
    vector<Eigen::Vector3d> collarVel;    
    vector<Eigen::Vector3d> collarPoints;    
    
    
    void print(lexer*);
    
    void buildNet_bag(lexer*);
	void buildNet_cyl(lexer*);
    void buildNet_wall(lexer*);
    
    void fillRhs_bag(lexer*);
    
    void fillRhs_Morison(lexer*);
    void morisonForceCoeff(double&, double&, const double&); 
 
    void fillRhs_Screen(lexer*);  
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
	double EA, w, rho_c, l_c, d_c;
	
	// Mesh
	double origin_x, origin_y, origin_z, phi, theta, psi;
	double *l0, *l;
    double L, D, al, ad, Fg, beta, gamma; 
    int nd, nl, niK, nbK, nK, nf;
    MatrixXd fi, A, B, Bh;
    double **fb, **K, **K_;
    int *Pb, *Nb, *Pi, *Ni;
    vector<vector<int> > meshID;
    
    // Net mesh
    
    vector<double> meshPoints;
    vector<int> meshSide;
    vector<double> tetVol;
    MatrixVd tri_x, tri_y, tri_z, tri_vel, tetC;
    MatrixVi tetMesh;
    
    // Raytracing
    fieldint5 cutl,cutr;
    double xs,xe,ys,ye,zs,ze;
    int tend;  
	
    // Reini
	vec f_,frk1,frk2,L_, dt;
	int reiniter;
	double xmin,xplus,ymin,yplus,zmin,zplus;
	double dstx,dsty,dstz,lsSig,dnorm,op,lsv,sign;

	// Forces
    double Fx,Fy,Fz;
	double **coupledField, **v_t, **v_n;
	int **nfK;

	// Print
	double *x, *y, *z, *T;
	char name[100];
	ofstream eTout;
	double printtime;
};

#endif
