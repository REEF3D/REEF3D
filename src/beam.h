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
                          
The solver was written in C++ by Blake Ashby (bmashby@stanford.edu) (version Nov 15, 2002). 
It is modified for C++ from the code RADAU5 originally written in FORTRAN (version July 9, 
1996, latest small correction: January 18, 2002) by E. Hairer (ernst.hairer@math.unige.ch) 
and G. Wanner (gerhard.wanner@math.unige.ch), Universite de Geneve. 

--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"vec.h"
#include"boundarycheck.h"
#include<iostream>
#include<vector>
#include <Eigen/Dense>

class lexer;

using namespace std;

#ifndef BEAM_H_
#define BEAM_H

class beam
{
public:
	
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    typedef Eigen::Matrix<double,3,Eigen::Dynamic> Matrix3Xd;
    typedef Eigen::Matrix<double,4,Eigen::Dynamic> Matrix4Xd;

    beam(int);
    ~beam();
    
    virtual void iniMaterial();
    virtual void meshBeam(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::Vector3d&);
    virtual void meshBeam(double, double, double, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&);
    virtual void setConstantLoads(Matrix3Xd&, Matrix4Xd&, const Matrix3Xd&, const Matrix3Xd&, const Matrix4Xd&, const Matrix4Xd&);
    virtual void setVariableLoads(Matrix3Xd&, Matrix4Xd&, const Matrix3Xd&, const Matrix3Xd&, const Matrix4Xd&, const Matrix4Xd&, const double);
    virtual void setFieldBC(Matrix3Xd&, Matrix3Xd&, Matrix4Xd&, Matrix4Xd&, Matrix4Xd&, Matrix3Xd&, Matrix4Xd&, Matrix3Xd&, double, int);
    virtual void print(lexer *p);

    void iniDamping(double, double, double, double, double, double, bool);
    void iniBeam(double, double, double, double, double, double, double, double, double);
    void iniSolver();
    void Integrate(double, double);

    void getTransPos(Matrix3Xd& c_){c_ = c;};
    void getTransVel(Matrix3Xd& cdot_){cdot_ = cdot;};
    void getRotPos(Matrix4Xd& q_){q_ = q;};
    void getRotVel(Matrix4Xd& qdot_){qdot_ = qdot;};
    
    Eigen::Vector3d getOmega(const Eigen::Vector4d&, const Eigen::Vector4d&);
    Eigen::Vector3d getOmega0(const Eigen::Vector4d&, const Eigen::Vector4d&);
    Eigen::Vector3d rotVec(const Eigen::Vector3d&, const Eigen::Vector4d&);

    double getTensLoc(int n){return f0.col(n).norm();};
    Eigen::Vector3d getTensGlob(int n){return R(q.col(n+1))*f0.col(n).tail(3);};
    
private:

    // Initialisation
    void iniFields();

    // Runtime
    void setState(double*);
    void setStateDot(double*);
    void getState(double*);
    void rhs(Matrix3Xd&, Matrix3Xd&, Matrix4Xd&, Matrix4Xd&, double);
    void resetSolver();

    // Inner forces and moments
    Eigen::Vector4d f0_(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector4d&, const Eigen::Vector4d&, const Eigen::Vector4d&);
    Eigen::Vector4d m0_(const Eigen::Vector4d&, const Eigen::Vector4d&, const Eigen::Vector4d&, const Eigen::Vector4d&, const Eigen::Vector4d&, const Eigen::Vector4d&);

    // Help routines
    Eigen::Vector4d qconj(const Eigen::Vector4d&);
    Eigen::Vector4d qMult(const Eigen::Vector4d&, const Eigen::Vector4d&);
    Eigen::Vector4d qMult(const Eigen::Vector4d&, const Eigen::Vector4d&, const Eigen::Vector4d&);
    Eigen::Vector4d qMult(int, const Eigen::Vector3d&, const Eigen::Vector4d&, const Eigen::Vector4d&);
    Eigen::Matrix3d R(const Eigen::Vector4d&);
    void calcQ(const Eigen::Vector4d&);
    void calcInvM(const Eigen::Vector4d&);

    // Solver
    int CoreIntegrator(double, double);
    double ContinuousOutput(unsigned);
    void ComputeJacobian();
    int DecompReal();
    int DecompComplex();
    int LinearSolve();
    int ErrorEstimate();
    void Function(double, double*, double*);
    void Jacobian(double, double*, double**);
    void Mass(double **M);
    int SolutionOutput();

        
    // ------ 
       
    // Beam number
    int nBeam;

    // Material
    double E, G, A, L, rho;
    Eigen::Matrix3d I, Ceps, Ckappa, Cepsdot, Ckappadot;
    Eigen::Matrix4d invIq, Iq, invM, Q, J;
    bool compression;

    // Mesh
    double dZ;
    int Ne;

    // Fields
    double *y;
    Matrix3Xd Fext, c, c0, cdot, cdotdot, rhs_cdot, f; 
    Matrix4Xd Mext, q, q0, qdot, rhs_qdot, f0, m0; 
    Eigen::Vector3d corr_k, mdot;
    Eigen::Vector4d qb, dcdz, dc0dz, dcdotdz, fdot, qmult, fq;
    Eigen::Vector4d dummy;

    // Print
	char name[100];
	double printtime;

    // Solver
    double **e1, **e2r, **e2i, **fjac, **fmas;
    double *rtoler, *atoler, *z1, *z2, *z3, *y0, *scal, *f1, *f2, *f3, *cont;
    double hinit, delta_t, dx, hmax, uround, safe, facl, facr, fnewt, quot1, quot2, thet, fac1, alphn, betan, err, hold, xold, x, x_end, xd;
    const int iout, imas, ijac;
    int *ip1, *ip2, *iphes;
    int n_dim, itoler, mljac, mujac, m1, m2, mlmas, mumas, nmax, nit, nind1, nind2, nind3, npred, nm1, ldjac, lde1, ldmas, ijob, njac, ndec, nsol, mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag, naccpt, nfcn, nstep, nrejct;
    bool startn, hess, pred, caljac, calhes, first, reject, implct, jband;

    // Solver routines
    template<class T>
    inline const T& max_(const T& a, const T& b)
    { return a > b ? a : b; }

    template<class T>
    inline const T& min_(const T& a, const T& b)
    { return a < b ? a : b; }

    // Matrix Triangularization by Gaussian Elimination
    int dec(const int n, double **A, int *ip);

    // Solution of linear system A*x = b
    void sol(const int n, double **A, double *b, int *ip);

    // Matrix Triangularization by Gaussian Elimination of a Hessenberg
    // matrix with lower bandwidth lb
    int dech(const int n, double **A, int lb, int *ip);

    // Solution of linear system A*x = b -- Hessenberg matrix
    void solh(const int n, double **A, int lb, double *b, int *ip);

    // Matrix Triangularization by Gaussian Elimination for complex matrices
    int decc(const int n, double **AR, double **AI, int *ip);

    // Solution of linear system A*x = b -- complex matrices
    void solc(const int n, double **AR, double **AI, double *br,
        double *bi, int *ip);

    // Matrix Triangularization by Gaussian Elimination -- Hessenberg, complex
    // matrices
    int dechc(const int n, double **AR, double **AI, int lb, int *ip);

    // Solution of linear system A*x = b -- Hessenberg, complex matrices
    void solhc(const int n, double **AR, double **AI, int lb,
        double *br, double *bi, int *ip);

    //Matrix Triangularization by Gaussian Elimination -- banded matrix
    int decb(const int n, double **A, int ml, int mu, int *ip);

    // Solution of linear system A*x = b -- banded matrix
    void solb(const int n, double **A, int ml, int mu, double *b, int *ip);

    //Matrix Triangularization by Gaussian Elimination -- banded, complex matrices
    int decbc(const int n, double **AR, double **AI, int ml, int mu, int *ip);

    // Solution of linear system A*x = b -- banded, complex matrices
    void solbc(const int n, double **AR, double **AI, int ml, int mu,
        double *br, double *bi, int *ip);

    // reduces a submatrix to upper Hessenberg form
    void elmhes(const int n, int low, int igh, double **A, int *inter);
};

#endif
