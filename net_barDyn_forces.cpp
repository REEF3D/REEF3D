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

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"
#include"vrans.h"


void net_barDyn::gravityForce(lexer *p)
{    
    for (int i = 0; i < nK; i++)
    {
        // Assign gravity force to knot
        forces_knot(i, 2) += -9.81*(mass_knot(i) - weight_knot(i));
    }
}


void net_barDyn::inertiaForce(lexer *p)
{   
    int index = 0;

    for (int i = 0; i < nK; i++)
    {
        if (i >= nfK[0][0])  // then i is an inner knot 
        {
            int*& barsiKI = nfK[index];
            int& kI = barsiKI[0];
                     
            // Assign inertia force to knot
            forces_knot(kI,0) += 2.0*weight_knot(kI)*(coupledField[kI][0] - coupledFieldn[kI][0])/dt_; 
            forces_knot(kI,1) += 2.0*weight_knot(kI)*(coupledField[kI][1] - coupledFieldn[kI][1])/dt_; 
            forces_knot(kI,2) += 2.0*weight_knot(kI)*(coupledField[kI][2] - coupledFieldn[kI][2])/dt_;   
    
            index++;
        }

    }
}


void net_barDyn::dragForce(lexer *p)
{    
    // Screen force model of Kristiansen (2012)
    // Assign hydrodynamic forces to knot from each adjoint screen
    
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    
    int nBars;
    int index = 0;
    
    Vector3d v_rel(1,3);
    Vector3d n_v(1,3);
    Vector3d Fadd(1,3);
    double v_mag, rho;
    
    for (int i = 0; i < nK; i++)
    {
        if (i >= nfK[0][0])  // then i is an inner knot 
        {
            int*& barsiKI = nfK[index];
            int& kI = barsiKI[0];
            
            //- Count number of screens
            
            nBars = 4;
            
            for (int k = 1; k < 5; k++)
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
        
            //- Add hydrodynamic force contributions from each screen
            
            if (nBars == 2)         // Corner screens
            {
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[2]);
            }      
            else if (nBars == 3)    // Edge screens
            {
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[3]);
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[2],barsiKI[3]);
            }
            else                    // Inner screens
            {
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[4]);
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[4],barsiKI[2]);
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[2],barsiKI[3]);
                forces_knot.row(i) += screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[3]);            
            }

            // Sinker drag
            if (i >= nK - nd - 1)
            {
                double Re = sinker_d*v_mag/p->W2;

                double logRe = log10(Re);

                double cd_circ = min(1.45 + 8.55*pow(Re, -0.9), 7.0);
                if (Re > 32)
                {
                    cd_circ =
                    -78.46675 + 254.73878*logRe - 327.8864*pow(logRe,2.0)
                    + 223.64577*pow(logRe,3.0) - 87.92234*pow(logRe,4.0)
                    + 20.00769*pow(logRe,5.0) - 2.44894*pow(logRe,6.0)
                    + 0.12479*pow(logRe,7.0);
                }
    
                Fadd = 0.5*rho*sinker_d*sinker_l*v_rel*v_mag*cd_circ;
 
                forces_knot.row(i) += Fadd;
 
                Fx += Fadd(0);
                Fy += Fadd(1);
                Fz += Fadd(2);
            }

            // Knot drag
            Fadd = 0.5*rho*PI/4.0*knot_d*knot_d*v_rel*v_mag*1.5;
            forces_knot.row(i) += Fadd;
            Fx += Fadd(0);
            Fy += Fadd(1);
            Fz += Fadd(2);

            index++;    
        }    
    }
}


Eigen::Vector3d net_barDyn::screenForce
(
    lexer *p,
    const double& rho,
    const Vector3d& v_rel,
    const Vector3d& n_d,
    const double& v_mag,
    const int b1,
    const int b2
)
{
    //- Normal vector of screen n_s
    
    Vector3d s1 = x_.row(Ni[b1]) - x_.row(Pi[b1]);
    
    Vector3d s2 = x_.row(Ni[b2]) - x_.row(Pi[b2]); 

    Vector3d n_s = s1.cross(s2);
    
    double A_panel = n_s.norm();

    n_s /= A_panel;
    n_s = SIGN(n_d.dot(n_s))*n_s;
    

    //- Area and solidity of screen (Lader, 2003)
    
    double As, Sn;

    if (b1 <= nd)   // add one screen to the upper row of meshes
    {
        As = 0.5*A_panel;

        Sn = p->X321_Sn[nNet]*0.5*l0[b1]*l0[b2]/As;
    }
    else
    {
        As = 0.25*A_panel;

        Sn = p->X321_Sn[nNet]*0.25*l0[b1]*l0[b2]/As;
    }


    //- Angle between velocity and normal vector
    
    double thetan = acos(n_d.dot(n_s));


    //- Normal vector of lift force
    
    Vector3d n_l = (n_d.cross(n_s)).cross(n_d).normalized();


    //- Get drag and lift force coefficients
    
    double cd, cl;

    double v_mag_corr = v_mag;
    
    double error = 1.0;
    int nIt = 0;

    while (error > 1e-3 && nIt < 10)
    {
        error = v_mag_corr;    
        
        screenForceCoeff(p,cd,cl,v_mag_corr,thetan,Sn);
        
        // Froude momentum theory
        v_mag_corr = v_mag*cd/(2.0*(sqrt(1.0 + cd) - 1.0)); 

        error = fabs(v_mag_corr - error);
        
        nIt++;
    }
    
    if (isnan(v_mag_corr))
    {
        v_mag_corr = v_mag;
        screenForceCoeff(p,cd,cl,v_mag_corr,thetan,Sn);
    }


    //- Calculate local forces
    
    double Fd = 0.5*rho*As*pow(v_mag_corr,2.0)*cd;
    double Fl = 0.5*rho*As*pow(v_mag_corr,2.0)*cl;
    
    
    //- Calculate global forces
    
    Fx += Fd*n_d(0) + Fl*n_l(0);
    Fy += Fd*n_d(1) + Fl*n_l(1);
    Fz += Fd*n_d(2) + Fl*n_l(2);

    return (Fd*n_d + Fl*n_l);
}


void net_barDyn::screenForceCoeff
(
    lexer *p,
    double& cd,
    double& cl,
    const double& v_mag,
    const double& thetan,
    const double& Sn
)
{
    //- Drag and lift force coefficients according to Kristiansen (2012)
    
    double a3 = 0.01387;    
    double a5 = 0.01361;  

    double b2 = 1.45;       //1.22905;  
    double b4 = 0.05369;    //0.11155;
    double b6 = 0.000367;   //0.0;

    double Re = d_c*v_mag/(p->W2*(1-Sn));

    double logRe = log10(Re);

    //double cd_circ = 1.1 + 4.0*pow(Re, -0.5);
    
    double cd_circ = min(1.45 + 8.55*pow(Re, -0.9), 7.0);
    if (Re > 32)
    {
        cd_circ =
        -78.46675 + 254.73878*logRe - 327.8864*pow(logRe,2.0)
        + 223.64577*pow(logRe,3.0) - 87.92234*pow(logRe,4.0)
        + 20.00769*pow(logRe,5.0) - 2.44894*pow(logRe,6.0)
        + 0.12479*pow(logRe,7.0);
    }
    

    double cd0 = cd_circ*Sn*(2.0 - Sn)/(2.0*pow(1.0 - Sn,2.0));

    double cn_piFourth =
        cd_circ*Sn*(2.0 - Sn)/(2.0*pow(1.0 - Sn,2.0))*pow(cos(PI/4.0),2.0);

    double ct_piFourth = PI/4.0*4.0*cn_piFourth/(8.0 + cn_piFourth);

    double cl0 = (0.5*cd0 - ct_piFourth)/sqrt(2.0);

    cd = cd0*((1.0 - a3 - a5)*cos(thetan) + a3*cos(3.0*thetan) + a5*cos(5.0*thetan));

    cl = cl0*(b2*sin(2.0*thetan) + b4*sin(4.0*thetan) + b6*sin(6.0*thetan)); 
}

void net_barDyn::netForces
(
    lexer *p,
	double& Xne, double& Yne, double& Zne,
	double& Kne, double& Mne, double& Nne
)
{
    Xne = 0.0;
	Yne = 0.0;
	Zne = 0.0;        
    Kne = 0.0;
	Mne = 0.0;
	Nne = 0.0;
    
    double Xne_i, Yne_i, Zne_i;
    Vector3d T_knot, x_ij;
    int knotJ, barI;

    Tne = 0.0;


    // Calculate force distribution from vertical twines at the top
    for (int knotI = 0; knotI < nbK; knotI++)
    {
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
        Xne_i = T_(barI)*x_ij(0);
        Yne_i = T_(barI)*x_ij(1);
        Zne_i = T_(barI)*x_ij(2);

        // Add force distribution
        Tne += T_(barI);
        Xne += Xne_i;
        Yne += Yne_i;
        Zne += Zne_i;

        // Add moment distribution
        Kne += (x_(knotI,1) - p->yg)*Zne_i - (x_(knotI,2) - p->zg)*Yne_i;
        Mne += (x_(knotI,2) - p->zg)*Xne_i - (x_(knotI,0) - p->xg)*Zne_i;
        Nne += (x_(knotI,0) - p->xg)*Yne_i - (x_(knotI,1) - p->yg)*Xne_i;
    } 
}
