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

#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"
#include"vrans.h"


void net_sheet::gravityForce(lexer *p)
{    
    for (int knotI = 0; knotI < nK; knotI++)
    {
        // Assign gravity force to knot
        forces_knot(knotI,2) -= 9.81*(mass_knot(knotI) - weight_knot(knotI));
    }
}


void net_sheet::inertiaForce(lexer *p)
{   
    if (dt_ < 1e-10) dt_ = 1e10;
    
    for (int knotI = 0; knotI < nK; knotI++)
    {
        // Assign inertia force to knot
        forces_knot(knotI,0) = 2.0*weight_knot(knotI)*(coupledField[knotI][0] - coupledFieldn[knotI][0])/dt_;
        forces_knot(knotI,1) = 2.0*weight_knot(knotI)*(coupledField[knotI][1] - coupledFieldn[knotI][1])/dt_;
        forces_knot(knotI,2) = 2.0*weight_knot(knotI)*(coupledField[knotI][2] - coupledFieldn[knotI][2])/dt_;   
    }
}


void net_sheet::dragForce(lexer *p)
{    
    // Assign hydrodynamic forces to knot 
    Vector3d side1, side2, normalVec;
    Vector3d v_rel, n_d, n_s, n_l;
    double mag, v_mag, rho, error, area, thetan, cd, cl, v_mag_corr;
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
    
    for (int knotI = 0; knotI < nK; knotI++)
    {
        //Calculate normal vector
		x0 = tri_x[knotI][0];
		x1 = tri_x[knotI][1];
		x2 = tri_x[knotI][2];
		
		y0 = tri_y[knotI][0];
		y1 = tri_y[knotI][1];
		y2 = tri_y[knotI][2];
		
		z0 = tri_z[knotI][0];
		z1 = tri_z[knotI][1];
		z2 = tri_z[knotI][2];  
        
        side1 << x1-x0, y1-y0, z1-z0;
        side2 << x2-x0, y2-y0, z2-z0;
        
        normalVec = side1.cross(side2);
        mag = normalVec.norm();
        area = 0.5*mag;
        normalVec /= mag;
        
        // Calculate relative velocity at knot
        v_rel << coupledField[knotI][0] - xdot_(knotI,0), 
                 coupledField[knotI][1] - xdot_(knotI,1), 
                 coupledField[knotI][2] - xdot_(knotI,2);
                     
        // Calculate normal velocity vector
        v_mag = v_rel.norm();
        n_d = v_rel/(v_mag + 1e-10);

        // Access density at knot
        rho = coupledField[knotI][3] > 900 ? coupledField[knotI][3] : 0.0;   
       
        // Correct direction of normal vector of triangle
        n_s = SIGN(n_d.dot(normalVec))*normalVec;
        
        // Angle between velocity and normal vector
        thetan = acos(n_d.dot(n_s));     

        // Normal vector of lift force
        n_l = (n_d.cross(n_s)).cross(n_d).normalized();
        
        //- Get drag and lift force coefficients
        v_mag_corr = v_mag;
        
        error = 1.0;
        int nIt = 0;

        while (error > 1e-3 && nIt < 10)
        {
            error = v_mag_corr;    
            
            screenForceCoeff(p,cd,cl,v_mag_corr,thetan,p->X321_Sn[nNet]);
            
            // Froude momentum theory
            v_mag_corr = v_mag*cd/(2.0*(sqrt(1.0 + cd) - 1.0)); 

            error = fabs(v_mag_corr - error);
            
            nIt++;
        }
        
        if (std::isnan(v_mag_corr))
        {
            v_mag_corr = v_mag;
            screenForceCoeff(p,cd,cl,v_mag_corr,thetan,p->X321_Sn[nNet]);
        }            
        

        // Save directional forces
        forces_knot.row(knotI) += 0.5*area*rho*pow(v_mag_corr,2.0)*(cd*n_d + cl*n_l);

        // Knot drag
        forces_knot.row(knotI) += 0.5*rho*PI/4.0*knot_d*knot_d*v_rel*v_mag*1.5;
    }
}


void net_sheet::screenForceCoeff
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
    

    if(p->Y1 == 1)
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

void net_sheet::netForces
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
    
    double FxI, FyI, FzI;
    
    for (int knotI = 0; knotI < nK; knotI++)
    {
        // Assign forces
        FxI = forces_knot(knotI, 0);
        FyI = forces_knot(knotI, 1);
        FzI = forces_knot(knotI, 2);

        // Add force distribution
        Xne += FxI; 
        Yne += FyI;
        Zne += FzI;

        // Add moment distribution
        Kne += (x_(knotI,1) - p->yg)*FyI - (x_(knotI,2) - p->zg)*FzI;
        Mne += (x_(knotI,2) - p->zg)*FzI - (x_(knotI,0) - p->xg)*FxI;
        Nne += (x_(knotI,0) - p->xg)*FxI - (x_(knotI,1) - p->yg)*FyI;
    } 
}
