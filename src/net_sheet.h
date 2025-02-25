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

#ifndef NET_SHEET_H_
#define NET_SHEET_H_

class net_sheet : public net, public boundarycheck
{
public:
    net_sheet(int, lexer*);
    virtual ~net_sheet();
    
    virtual void start(lexer*, fdm*, ghostcell*,double,Eigen::Matrix3d);
    virtual void initialize(lexer*, fdm*, ghostcell*);
    virtual void netForces(lexer*, double&, double&, double&, double&, double&, double&);
    virtual const EigenMat& getLagrangePoints(){return lagrangePoints;} 
    virtual const EigenMat& getLagrangeForces(){return lagrangeForces;} 
    virtual const EigenMat& getCollarVel(){return collarVel;} 
    virtual const EigenMat& getCollarPoints(){return collarPoints;} 

    
private:

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorXd;
    typedef Eigen::Matrix<double, 3, 3> Matrix3d;
    typedef Eigen::Matrix<double, 1, 3> Vector3d;
    
    typedef vector<vector<double> > MatrixVd;
    typedef vector<vector<int> > MatrixVi;
    
    
    // Preprocessing
    void ini(lexer*, fdm*, ghostcell*); 
    void rotation_tri(lexer*,double,double,double,double&,double&,double&, const double&, const double&, const double&);
    
    // Runtime
    void updateField(lexer*, fdm*, ghostcell*, int);
    void print(lexer*);
    Eigen::VectorXd timeWeight(lexer*);
    void gravityForce(lexer*);
    void dragForce(lexer*);
    void inertiaForce(lexer*);
    void screenForceCoeff(lexer*,double&, double&, const double&, const double&, const double&);
   


    void vransCoupling(lexer*, fdm*, ghostcell*);
    void triangulation(lexer*, fdm*, ghostcell*);
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
    double rho_c, l_c, d_c;
    
    // Mesh
    int nK;
    MatrixXd x0_, x_, xdot_;
    double dt_; 
    VectorXd mass_knot, weight_knot, added_mass;
    MatrixXd forces_knot;
    double **coupledField, **coupledFieldn;
    
    // Net mesh
    int tend;  
    MatrixVd tri_x, tri_y, tri_z;
    MatrixVd tri_x0, tri_y0, tri_z0;
    vector<Eigen::Vector3d> lagrangePoints;    
    vector<Eigen::Vector3d> lagrangeForces;  
    vector<Eigen::Vector3d> collarVel;    
    vector<Eigen::Vector3d> collarPoints;    

    // Forces
    double Fx,Fy,Fz;

    // Knots
    double knot_d;

    // Probe Points
    Eigen::VectorXi probeKnot;

    // Print
    char name[100];
    double printtime;
};

#endif
