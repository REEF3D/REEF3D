/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"ship_models.h"
#include<cmath>

static const double PI_SHIP = 3.14159265358979323846;

const double ship_models::Re_min = 1.0e5;

double ship_models::cf_ittc57(double Re)
{
    const double R = Re>Re_min ? Re : Re_min;
    const double lg = log10(R) - 2.0;
    
    return 0.075/(lg*lg);
}

double ship_models::friction(double rho, double nu, double S, double L, double k, double u, double &Re, double &CF)
{
    Re = fabs(u)*L/nu;
    CF = cf_ittc57(Re);
    
    return -0.5*rho*S*(1.0+k)*CF*fabs(u)*u;
}

void ship_models::crossflow(double rho, double Cd, const std::vector<double> &xs, const std::vector<double> &dx, const std::vector<double> &T,
                            double v, double r, double &Y, double &N)
{
    Y = N = 0.0;
    
    for(size_t s=0; s<xs.size(); ++s)
    {
        const double vl = v + xs[s]*r;   // local transverse velocity of the strip
        const double dY = -0.5*rho*Cd*T[s]*fabs(vl)*vl*dx[s];
        
        Y += dY;
        N += xs[s]*dY;
    }
}

double ship_models::roll_damping(double B44, double B44q, double p)
{
    return -B44*p - B44q*fabs(p)*p;
}

void ship_models::propeller(double rho, double n, double D, const double *kt, const double *kq, double Va,
                            double &J, double &KT, double &KQ, double &T, double &Q)
{
    J = fabs(n)>1.0e-12 ? Va/(n*D) : 0.0;
    
    KT = kt[0] + kt[1]*J + kt[2]*J*J;
    KQ = kq[0] + kq[1]*J + kq[2]*J*J;
    
    T = rho*n*fabs(n)*pow(D,4.0)*KT;
    Q = rho*n*fabs(n)*pow(D,5.0)*KQ;
}

void ship_models::rudder_mmg(double rho, const rudder_param &R, double u, double v, double r, double delta,
                             double D, double n, double KT, double uP,
                             double &X, double &Y, double &N, double &K, double &alphaR, double &UR, double &FN)
{
    const double U = sqrt(u*u + v*v);
    const double HR = sqrt(R.AR*R.Lambda);
    const double falpha = R.falpha>0.0 ? R.falpha : 6.13*R.Lambda/(R.Lambda + 2.25);
    
    // longitudinal inflow: propeller slipstream on the fraction eta = D/HR of the rudder span,
    // written with KT n^2 D^2 = KT J^2 uP^2 / ... so that it holds at J = 0 (bollard)
    double uR;
    
    if(D>0.0)
    {
        const double eta = D/HR<1.0 ? D/HR : 1.0;
        const double c = uP*uP + 8.0*KT*n*n*D*D/PI_SHIP;
        const double us = uP + R.kappa*(sqrt(c>0.0 ? c : 0.0) - uP);
        uR = R.eps*sqrt(eta*us*us + (1.0 - eta)*uP*uP);
    }
    else
    uR = R.eps*uP;
    
    // lateral inflow with flow straightening; MMG frame (y to starboard): v_M = -v, r_M = -r,
    // beta_M = atan(-v_M/u), beta_R = beta_M - lR r_M/U
    double vR = 0.0;
    
    if(U>1.0e-10)
    {
        const double betaR = atan2(v,u) + R.lR*r/U;
        const double gR = (R.gammaRp>=0.0 && betaR>0.0) ? R.gammaRp : R.gammaR;
        vR = U*gR*betaR;
    }
    
    alphaR = delta - atan2(vR, uR>1.0e-10 ? uR : 1.0e-10);
    
    UR = sqrt(uR*uR + vR*vR);
    FN = 0.5*rho*R.AR*UR*UR*falpha*sin(alphaR);
    
    // MMG forces, Y and N with the sign of the ship frame (y to port, z up)
    X = -(1.0 - R.tR)*FN*sin(delta);
    Y =  (1.0 + R.aH)*FN*cos(delta);
    N =  (R.xR + R.aH*R.xH)*FN*cos(delta);
    K = -R.zR*Y;
}


void ship_models::mmg_hull(double rho, double L, double d, const double *c, double u, double vm, double r, double Umin,
                           double &X, double &Y, double &N)
{
    // forces written as polynomials in (u, vm, r) with U^2 v' = U vm, U^2 r' = U r L, ... so that
    // they stay finite at small speed
    const double U = sqrt(u*u + vm*vm);
    const double Ud = U>Umin ? U : Umin;
    const double rl = r*L;
    
    const double fX = -c[0]*U*U + c[1]*vm*vm + c[2]*vm*rl + c[3]*rl*rl + c[4]*vm*vm*vm*vm/(Ud*Ud);
    const double fY = U*(c[5]*vm + c[6]*rl) + (c[7]*vm*vm*vm + c[8]*vm*vm*rl + c[9]*vm*rl*rl + c[10]*rl*rl*rl)/Ud;
    const double fN = U*(c[11]*vm + c[12]*rl) + (c[13]*vm*vm*vm + c[14]*vm*vm*rl + c[15]*vm*rl*rl + c[16]*rl*rl*rl)/Ud;
    
    X = 0.5*rho*L*d*fX;
    Y = 0.5*rho*L*d*fY;
    N = 0.5*rho*L*L*d*fN;
}

double ship_models::mmg_wake(double wP0, double C1, double C2p, double C2n, double betaP)
{
    const double C2 = betaP>0.0 ? C2p : C2n;
    
    return 1.0 - (1.0 - wP0)*(1.0 + (1.0 - exp(-C1*fabs(betaP)))*(C2 - 1.0));
}
