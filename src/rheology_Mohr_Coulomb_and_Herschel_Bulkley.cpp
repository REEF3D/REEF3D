/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

// Three phase viscosity - air, slurry and gravel
// Air phase is covered in fluid_update_rheology.cpp
// inputs: n=W98, tau_00 & C* & b & m_y=W106
double rheology_f::Mohr_Coulomb_and_Herschel_Bulkley(lexer* p, fdm* a, ghostcell* pgc)
{
    // fines + suspension = slurry: Herschel–Bulkley
    const double shear_rate = strainterm(p,a);
    const double b = p->W106_b;
    const double n = p->W98; // calibration : shear thinning vs shear thickening
    const double C = p->W105_C_total; // volumetric solid concentration: the volume of all solid particles including fine material relative to the volume of the debris-flow material including water m^3/m^3
    const double tau_00 = p->W105_tau_00; // calibration : grid resolution sensitivity of the shear rate
    const double tau_0 = C<0.47?tau_00:tau_00*exp(5*(C-0.47));
    const double C_kaolinite_chlorite = p->W105_C_kaolinite_chlorite;
    const double C_illite = p->W105_C_illite;
    const double C_montmorillonite = p->W105_C_montmorillonite;
    const double P_0 = C_kaolinite_chlorite + 1.3*C_illite + 1.7*C_montmorillonite;
    const double P_1 = P_0>0.27 ? 0.7*P_0 : P_0;
    const double tau_y = tau_0*C*C*exp(22*C*P_1);
    const double k = b*tau_y; // d_max<0.4mm
    double mu_2 = k * pow(fabs(shear_rate),n-1) + tau_y*pow(fabs(shear_rate),-1);
    mu_2 = MIN(mu_2,p->W95);

    // gravel: Coulomb viscoplastic
    const double mu_min = p->W107_mu_min; // minimal dynamic viscosity
    const double delta = (p->W107_delta/180.0)*M_PI; // internal friction angle
    const double m_y = p->W106_m_y; // ]0,1]
    double mu_3 = mu_min + (shear_rate!=0?a->press(i,j,k)*sin(delta)/shear_rate*(1-exp(-m_y*shear_rate)):0.0);
    mu_3 = MIN(mu_3,p->W107_mu_0); // ToDo: check limiter

    const double a_2 = p->W108_a_2;
    const double a_3 = 1.0-a_2;
    // W1 is assumed to be the mean density of the mixture
    return (a_2*mu_2+a_3*mu_3)/p->W1;
}