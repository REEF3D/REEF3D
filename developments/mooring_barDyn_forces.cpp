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
--------------------------------------------------------------------*/

#include"mooring_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"
#include"vrans.h"


void mooring_barDyn::gravityForce(lexer *p)
{    
    for (int i = 0; i < nK; i++)
    {
        // Assign gravity force to knot
        forces_knot(i, 2) += -9.81*(mass_knot(i) - weight_knot(i));
    }
}


void mooring_barDyn::inertiaForce(lexer *p)
{   
    for (int index = 1; index < nK-1; index++)
    {
        int*& barsiKI = nfK[index];
        int& kI = barsiKI[0];
                 
        // Assign inertia force to knot
        forces_knot(kI,0) += 2.0*weight_knot(kI)*(coupledField[kI][0] - coupledFieldn[kI][0])/p->dt; 
        forces_knot(kI,1) += 2.0*weight_knot(kI)*(coupledField[kI][1] - coupledFieldn[kI][1])/p->dt; 
        forces_knot(kI,2) += 2.0*weight_knot(kI)*(coupledField[kI][2] - coupledFieldn[kI][2])/p->dt;   
    }
}


void mooring_barDyn::dragForce(lexer *p)
{    
    // Assign hydrodynamic forces to knot from each adjoint bar
    
    int nBars;
    
    Vector3d v_rel(1,3);
    Vector3d n_v(1,3);
    double v_mag, rho;
    
    for (int i = 1; i < nK-1; i++)
    {
        int*& barsiKI = nfK[i];
        int& kI = barsiKI[0];
        
        //- Count number of bars 
        
        nBars = 2;
        
        for (int k = 1; k < 3; k++)
        {
            if (barsiKI[k] == -1) nBars--;
        }
        
        
        //- Calculate relative velocity at knot
        
        v_rel << coupledField[kI][0] - xdot_(kI,0), 
                 coupledField[kI][1] - xdot_(kI,1), 
                 coupledField[kI][2] - xdot_(kI,2);
                 
        v_mag = v_rel.norm();
        n_v = v_rel/(v_mag + 1e-10);


        //- Access density at knot
        
        rho = coupledField[kI][3];   

        //- Add hydrodynamic force contributions from each bar
        
        double cd_t = 0.5;
        double cd_n = 2.5;
        double cm_t = 0.0;
        double cm_n = 3.8;

        double vx = v_rel(0); 
        double vy = v_rel(1); 
        double vz = v_rel(2);

        Vector3d x_ij = (x_.row(i) - x_.row(i-1)).normalized();
        
        double dl = (x_.row(i) - x_.row(i-1)).norm();

        double vt = vx*x_ij(0) + vy*x_ij(1) + vz*x_ij(2);

        double vnx = vx - vt*x_ij(0);
        double vny = vy - vt*x_ij(1);
        double vnz = vz - vt*x_ij(2);
        double vn_mag = sqrt(vnx*vnx + vny*vny + vnz*vnz);

        forces_knot(i,0) += 0.5*rho*d_c*dl/2.0*(cd_t*fabs(vt)*vt*x_ij(0) + cd_n*vn_mag*vnx);
        forces_knot(i,1) += 0.5*rho*d_c*dl/2.0*(cd_t*fabs(vt)*vt*x_ij(1) + cd_n*vn_mag*vny);
        forces_knot(i,2) += 0.5*rho*d_c*dl/2.0*(cd_t*fabs(vt)*vt*x_ij(2) + cd_n*vn_mag*vnz);
        
        x_ij = (x_.row(i+1) - x_.row(i)).normalized();
        
        dl = (x_.row(i+1) - x_.row(i)).norm();

        vt = vx*x_ij(0) + vy*x_ij(1) + vz*x_ij(2);

        vnx = vx - vt*x_ij(0);
        vny = vy - vt*x_ij(1);
        vnz = vz - vt*x_ij(2);
        vn_mag = sqrt(vnx*vnx + vny*vny + vnz*vnz);

        forces_knot(i,0) += 0.5*rho*d_c*dl/2.0*(cd_t*fabs(vt)*vt*x_ij(0) + cd_n*vn_mag*vnx);
        forces_knot(i,1) += 0.5*rho*d_c*dl/2.0*(cd_t*fabs(vt)*vt*x_ij(1) + cd_n*vn_mag*vny);
        forces_knot(i,2) += 0.5*rho*d_c*dl/2.0*(cd_t*fabs(vt)*vt*x_ij(2) + cd_n*vn_mag*vnz);
    }
}


void mooring_barDyn::bottomForce(lexer *p)
{
	double Kg = 3.0e9;
	double xi = 1.0;
    double dz;

    for (int i = 0; i < nK; i++)
    {
        dz = p->X311_zs[line] - x_(i,2);
        
        if (dz >= 0.0)
        {
            forces_knot(i, 2) += (Kg*d_c*dz - 2.0*xi*sqrt(Kg*mass_knot(i)*d_c)*max(xdot_(i,2),0.0));
        }
    }
}


void mooring_barDyn::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
    Vector3d T_knot, x_ij;
    int knotJ, barI;

    // Calculate force distribution from upper element
    int knotI = nK - 1;

    barI = nfbK[knotI];
    // Find bar vector
    if (Pi[barI] == knotI)
    {
        knotJ = Ni[barI];
    }
    else
    {
        knotJ = Pi[barI];
    }
    
    // Normal vector
    x_ij = (x_.row(knotJ) - x_.row(knotI)).normalized();

    // Calculate force contributions
    Xme = T_(barI)*x_ij(0);
    Yme = T_(barI)*x_ij(1);
    Zme = T_(barI)*x_ij(2);
}
