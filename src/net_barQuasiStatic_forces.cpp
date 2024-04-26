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
--------------------------------------------------------------------*/

#include"net_barQuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"
#include"vrans.h"


Eigen::Vector3d net_barQuasiStatic::screenForce
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
    // Normal vector of screen n_s
    Vector3d s1 = fi.row(b1);
    Vector3d s2 = fi.row(b2);

    Vector3d n_s = s1.cross(s2).normalized();

    n_s = SIGN(n_d.dot(n_s))*n_s;

    // Area and solidity of screen (Lader, 2003)
    double As, Sn;

    if (b1 <= nd)   // add one screen to the upper row of meshes
    {
        As = 0.5*l[b1]*l[b2]*s1.cross(s2).norm();

        Sn = p->X321_Sn[nNet]*0.5*l0[b1]*l0[b2]/As;
    }
    else
    {
        As = 0.25*l[b1]*l[b2]*s1.cross(s2).norm();

        Sn = p->X321_Sn[nNet]*0.25*l0[b1]*l0[b2]/As;
    }

    // Angle between velocity and normal vector
    double thetan = acos(n_d.dot(n_s));

    // Normal vector of lift force
    Vector3d n_l = (n_d.cross(n_s)).cross(n_d).normalized();

    // Get drag and lift force coefficients
    double cd, cl;

    double v_mag_corr = v_mag;
    
    double error = 1.0;
    int nIt = 0;

    while (error > 1e-3 && nIt < 10)
    {
        error = v_mag_corr;    
     
        screenForceCoeff(p,cd,cl,v_mag_corr,thetan,Sn);
        
       // v_mag_corr = v_mag/(1.0 - 0.25*cd);
        v_mag_corr = v_mag*cd/(2*(-1 + sqrt(1.0 + cd))); // Momentum theory

        error = fabs(v_mag_corr - error);
        
        nIt++;
    }
    
    if (std::isnan(v_mag_corr))
    {
        v_mag_corr = v_mag;
        screenForceCoeff(p,cd,cl,v_mag_corr,thetan,Sn);
    }

    // Calculate local forces
    double Fd = 0.5*rho*As*pow(v_mag_corr,2.0)*cd;
    double Fl = 0.5*rho*As*pow(v_mag_corr,2.0)*cl;
    
    // Calculate global forces
    Fx += Fd*n_d(0) + Fl*n_l(0);
    Fy += Fd*n_d(1) + Fl*n_l(1);
    Fz += Fd*n_d(2) + Fl*n_l(2);

    return (Fd*n_d + Fl*n_l);
}


void net_barQuasiStatic::screenForceCoeff
(
    lexer *p,
    double& cd,
    double& cl,
    const double& v_mag,
    const double& thetan,
    const double& Sn
)
{
    // Drag and lift force coefficients from Kristiansen (2012)
    
    double a3 = 0.0129;
    double a5 = 0.0151; 
    double a7 = 0.0;

    double b2 = 1.2044; 
    double b4 = 0.1466;
    double b6 = 0.135; 

    double Re = d_c*v_mag/(p->W2*(1-Sn));

    double logRe = log10(Re);

    double cd_circ =
        -78.46675 + 254.73878*logRe - 327.8864*pow(logRe,2.0)
        + 223.64577*pow(logRe,3.0) - 87.92234*pow(logRe,4.0)
        + 20.00769*pow(logRe,5.0) - 2.44894*pow(logRe,6.0)
        + 0.12479*pow(logRe,7.0);
   
    cd_circ = 1.1 + 4.0*pow(Re,(-0.5));
    

    double cd0 = cd_circ*Sn*(2.0 - Sn)/(2.0*pow(1.0 - Sn,2.0));

    double cn_piFourth =
        cd_circ*Sn*(2.0 - Sn)/(2.0*pow(1.0 - Sn,2.0))*pow(cos(PI/4.0),2.0);

    double ct_piFourth = PI/4.0*4.0*cn_piFourth/(8.0 + cn_piFourth);

    double cl0 = (0.5*cd0 - ct_piFourth)/sqrt(2.0);

    cd = cd0*((1.0 - a3 - a5)*cos(thetan) + a3*cos(3.0*thetan) + a5*cos(5.0*thetan));

    cl = cl0*(b2*sin(2.0*thetan) + b4*sin(4.0*thetan) + b6*sin(6.0*thetan));


    if(p->Y1==1)
    {
        // Simulation-based screen force model 

        double theta = thetan*180/PI;
        vector<double> p0 {-0.132, 340.797, -59.643, -9.129, 2.245, -12473.957, 0.063, 27831.591, 1458.245, 0.619};
        vector<double> p45 {-0.123,205.486,-40.789,-22.592,-10.21,-12828.831,0.297,60436.89,1787.102,-0.121,-24020.784,706.455,3775.06,253.464,-9.623,128457.65,-0.025,-21641.925,-9227634.55,0.368};
        vector<double> pl {-0.063,66.287,-10.840,-11.374,1.386,-2605.979,0.036,9838.141,267.245,0.113};

        cd0 = p0[0]*v_mag + p0[1]*d_c + p0[2]*l_c + p0[3]*v_mag*d_c + p0[4]*v_mag*l_c + p0[5]*d_c*l_c + p0[6]*v_mag*v_mag + p0[7]*d_c*d_c + p0[8]*l_c*l_c + p0[9];
    
        double cd45 = p45[0]*v_mag+p45[1]*d_c+p45[2]*l_c+p45[3]*v_mag*d_c+p45[4]*v_mag*l_c+p45[5]*d_c*l_c+p45[6]*v_mag*v_mag+p45[7]*d_c*d_c+p45[8]*l_c*l_c
                      + p45[9]*v_mag*v_mag*v_mag+p45[10]*l_c*l_c*l_c+p45[11]*v_mag*d_c*l_c
            +p45[12]*v_mag*d_c*d_c+p45[13]*v_mag*l_c*l_c+p45[14]*d_c*v_mag*v_mag+p45[15]*d_c*l_c*l_c+p45[16]*l_c*v_mag*v_mag+p45[17]*l_c*d_c*d_c+p45[18]*d_c*d_c*d_c+p45[19];	
    
        double cl45 = pl[0]*v_mag + pl[1]*d_c + pl[2]*l_c + pl[3]*v_mag*d_c + pl[4]*v_mag*l_c + pl[5]*d_c*l_c + pl[6]*v_mag*v_mag + pl[7]*d_c*d_c + pl[8]*l_c*l_c + pl[9];
    
        
        double a = -0.000493827160494*(cd45 - cd0) - 0.000246913580247*cd0; 
        double b = 0.044444444444444*(cd45 - cd0) + 0.011111111111111*cd0;
        cd = cd0 + a*theta*theta + b*theta;

        a = -0.000493827160494*cl45; 
        b = 0.044444444444444*cl45;
        cl = a*theta*theta + b*theta;
    }
}


void net_barQuasiStatic::morisonForceCoeff
(
    double& cn,
    double& ct,
    const double& vn_mag
)
{
    double Re, s;

    Re = d_c*vn_mag/(1.0e-6);

    // Tangential coefficient according to Tsukrov (2003)
    ct = PI*1.05e-3*(0.55*sqrt(Re) + 0.084*pow(Re,2.0/3.0));

    // Normal coefficient according to Tsukrov (2003)
    if (Re <= 1.0)
    {
        s = -0.077215665 + log(8.0/Re);

        cn = 8.0*PI/(Re*s)*(1.0 - 0.87*pow(s,-2.0));
    }
    else if (Re <= 30.0)
    {
        cn = 1.45 + 8.55*pow(Re,-0.9);
    }
    else
    {
        cn = 1.1 + 4.0*pow(Re,-0.5);
    }
}

void net_barQuasiStatic::netForces
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
}

