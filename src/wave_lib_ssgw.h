/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: 
    Csaba Pakozdi, Sébastien Fouques: Routine implementation
    Tobias Martin: Interface implementation

Based on Clamond and Dutykh (2018). Accurate fast computation of steady two-dimensional 
surface gravity waves in arbitrary depth. Journal of Fluid Mechanics, Vol. 844, pp. 491-518.
--------------------------------------------------------------------*/

#ifndef WAVE_LIB_SSGW_H_
#define WAVE_LIB_SSGW_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"increment.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <iostream>
#include <limits>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>

using namespace std;

struct Parameters
{
    int N, iter;
    double tol, res;
    
    double phaseVelocity, celerityCs, groupVelocity, BernoulliConstant;
    double waveHeight, creastHeight, troughHeight, waveLength, waveNumber, waterDepth, gravityAcc;
    double Impulse, PotentialEnergy, KineticEnergy, RadiationStress;
    double MomentumFlux, EnergyFlux;
};

class wave_lib_ssgw : public wave_lib_precalc, public wave_lib_parameters, public increment
{
public:
    wave_lib_ssgw(lexer*, ghostcell*);
	virtual ~wave_lib_ssgw();
    
    double wave_horzvel(lexer*,double,double,double);
    
    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_fi(lexer*,double,double,double);
    
    virtual void parameters(lexer*,ghostcell*);
    virtual void wave_prestep(lexer*,ghostcell*);
    
private:
    double singamma,cosgamma;
    
    void setWave(double k, double d, double H);
    void getPhysicsParameters();
    bool computeSurfaceVariables();
    void writeResult(const std::string folderName);
    void computeVelocityField(std::vector<double>& x, std::vector<double>& y, std::vector<double>& u, std::vector<double>& v);
    void computePotentialField(std::vector<double>& x, std::vector<double>& y, std::vector<double>& phi);
    Parameters ParameterValue;
    std::vector<double> xs, ys, us, vs, phis;
    
    bool resizing();
    void defineInput2Iteration(double k, double d, double H);
    void PetviashviliIteration();
    double modulo(double, double);
    
    int N, iter;
    double tol, err;
    bool allocated, infiniteDepth, surfaceCalculated;
    double H, k, d, g, L, lambda, c02, dalpha, dk;
    double B, ce, cs, cg, a, b, errfun;
    double intImpulse, intPotentialEnergy, intKineticEnergy;
    double intRadiationStress, intMomentumFlux, intEnergyFlux;
    double const GRAV = 9.81;
    double scaleVel, scaleLength, kd, kH2;
    Eigen::FFT<double> fft;

    double xcurr, xl, xr, yl, yr, fi, eta;
    
    Eigen::VectorXd dummy_real;
    Eigen::VectorXd Upsilon, CYs, C_hat;
    Eigen::VectorXd S2_hat, Upsilon2,CUpsilon;
    Eigen::VectorXd CUpsilon2,Cinf_hat,CIC_hat;
    Eigen::VectorXd L_hat,IL_hat,LUps,Nupsilon;
    Eigen::VectorXd IC, ydx, ICydx;
    Eigen::VectorXd U, Xs, Ys;
    Eigen::VectorXd vectorAlpha, vectorK;
    
    Eigen::VectorXcd zs, dzs, ws;
    Eigen::VectorXcd Zs, dZs, dummy_hat, dummy;
    Eigen::VectorXcd Ys_hat, Upsilon_hat,Upsilon2_hat;
    Eigen::VectorXcd CUps_hat,Nupsilon_hat,IH_hat;
};

#endif
