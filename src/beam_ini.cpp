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
--------------------------------------------------------------------*/

#include"beam.h"
#include<sys/stat.h>

void beam::iniBeam(double Ne_, double E_, double A_, double rho_, double L_, double G_, double IX_, double IY_, double IZ_)
{
    Ne = Ne_; 
    E = E_; 
    A = A_; 
    rho = rho_; 
    L = L_; 
    G = G_;
    
    // Set inertia matrix
    I << IX_, 0.0, 0.0, 
         0.0, IY_, 0.0, 
         0.0, 0.0, IZ_;
    
    // Initialise fields
    iniFields();

    // Initialise printing
    mkdir("./REEF3D_CFD_Beam",0777);
	printtime = 0.0;
}
    

void beam::iniMaterial()
{
    corr_k << 1.0, 1.0, 1.0;
    
    // Linear visco-elasticity model
    Ceps << E*A, 0, 0, 
            0, corr_k(0)*G*A, 0,
            0, 0, corr_k(1)*G*A;

    Ckappa << corr_k(2)*G, 0, 0,
              0, E, 0,
              0, 0, E;
    Ckappa *=I; 

    Iq(1,1) = I(0,0); Iq(2,2) = I(1,1); Iq(3,3) = I(2,2);
    invIq(1,1) = 1.0/I(0,0); invIq(2,2) = 1.0/I(1,1); invIq(3,3) = 1.0/I(2,2);
}

void beam::iniDamping(double cdx, double cdy, double cdz, double ckx, double cky, double ckz, bool compression_)
{
    // Linear visco-elasticity model
    Cepsdot << cdx, 0, 0,
               0, cdy, 0,
               0, 0, cdz; 

    Ckappadot << ckx, 0, 0,
                 0, cky, 0,
                 0, 0, ckz;

    // Compression effects on/off
    compression = compression_;
}


void beam::iniFields()
{
    // State vector
	n_dim = 6*(Ne+1) + 8*(Ne+2);
	y = new double[n_dim];
    
    // Matrices
    Iq = Eigen::Matrix4d::Zero(4,4);
    invIq = Eigen::Matrix4d::Zero(4,4);
    invM = Eigen::Matrix4d::Zero(4,4);

    // Edge matrices
    Fext = Matrix3Xd::Zero(3,Ne+1);
    c  = Matrix3Xd::Zero(3,Ne+1); 
    c0 = Matrix3Xd::Zero(3,Ne+1); 
    cdot  = Matrix3Xd::Zero(3,Ne+1); 
    cdotdot  = Matrix3Xd::Zero(3,Ne+1); 
    rhs_cdot = Matrix3Xd::Zero(3,Ne+1); 

    // Centre matrices with ghostpoints
    Mext = Matrix4Xd::Zero(4,Ne+2);
    q  = Matrix4Xd::Zero(4,Ne+2); 
    q0 = Matrix4Xd::Zero(4,Ne+2); 
    qdot  = Matrix4Xd::Zero(4,Ne+2); 
    rhs_qdot = Matrix4Xd::Zero(4,Ne+2);

    // Internal forces and moment matrices
    f  = Matrix3Xd::Zero(3,Ne+2);
    f0 = Matrix4Xd::Zero(4,Ne+2);
    m0 = Matrix4Xd::Zero(4,Ne+2);

    // Vectors
    dcdz = Eigen::Vector4d::Zero(4);
    dc0dz = Eigen::Vector4d::Zero(4);
    dcdotdz = Eigen::Vector4d::Zero(4);

    dummy = Eigen::Vector4d::Zero(4);

    // Bool
    compression = false;
}

void beam::meshBeam(double x_ini, double y_ini, double z_ini, Eigen::Vector3d& d1, Eigen::Vector3d& d2, Eigen::Vector3d& d3)
{
    // Direction vectors
    d1.normalize();
    d2.normalize();
    d3.normalize();
            
    double sp = (d3(1)*d2(2)-d3(2)*d2(1))*d1(0) +(d3(2)*d2(0)-d3(0)*d2(2))*d1(1)+(d3(0)*d2(1)-d3(1)*d2(0))*d1(2);
    if (sp < 0.0) d3 = -d3; 
    
    if( (1+d3(0)+d2(1)+d1(2)) < 1e-3)
    {
        d3 = -d3;
        d2 = -d2;
    }
    
    // Coordinates
    dZ = L/Ne; 
    for (int i=0; i < Ne+1; i++)
    {
        c.col(i) << x_ini + d1(0), y_ini + d1(1), z_ini + dZ*i*d1(2);
    }
    c0 = c; 

    // Quaternions
    for (int i=0; i < Ne+2; i++)
    {
        q(0,i) = sqrt( 1 + d3(0) + d2(1) + d1(2) )/2;
        q(1,i) = ( d1(1) - d2(2) )/(4*q(0,i));
        q(2,i) = ( d3(2) - d1(0) )/(4*q(0,i));
        q(3,i) = ( d2(0) - d3(1) )/(4*q(0,i));
        q.col(i)= qconj(q.col(i));
    }
    q0=q;

    // Fixed end quaternion
    qb = q.col(0);
}

void beam::meshBeam(const Eigen::VectorXd& x_, const Eigen::VectorXd& y_, const Eigen::VectorXd& z_, const Eigen::Vector3d& d0)
{
    // Coordinates
    c.row(0) = x_;
    c.row(1) = y_;
    c.row(2) = z_;
    c0 = c; 

    // Quaternions
    // R*q rotates point into body-fixed coordinate system
    for (int i=1; i < Ne+1; i++)
    {
        Eigen::Vector3d dc = (c.col(i) - c.col(i-1)).normalized(); 
        Eigen::Vector3d v = d0.cross(dc);
        double w = 1.0 + d0.dot(dc);
        q.col(i) << w, v(0), v(1), v(2);
        q.col(i).normalize();
    }
    q.col(0) = q.col(1); q.col(Ne+1) = q.col(Ne);
    q0=q;

    dZ = L/Ne; 

    // Fixed end quaternion
    qb << 1.0, 0.0, 0.0, 0.0;
}

void beam::iniSolver()
{
    setState(y);
    resetSolver();
	
    z1 = new double[n_dim];
	z2 = new double[n_dim];
	z3 = new double[n_dim];
	y0 = new double[n_dim];
	scal = new double[n_dim];
	f1 = new double[n_dim];
	f2 = new double[n_dim];
	f3 = new double[n_dim];
	cont = new double[4*n_dim];
	ip1 = new int[nm1];
	ip2 = new int[nm1];
	iphes = new int[n_dim];

	fjac = new double*[ldjac];
	fmas = new double*[ldmas];
	e1 = new double*[lde1];
	e2r = new double*[lde1];
	e2i = new double*[lde1];

	for (int i = 0; i < ldjac; i++) fjac[i] = new double[n_dim];

	for (int i = 0; i < ldmas; i++) fmas[i] = new double[n_dim];

	for (int i = 0; i < lde1; i++) {
		e1[i] = new double[nm1];
		e2r[i] = new double[nm1];
		e2i[i] = new double[nm1];
	}
}
