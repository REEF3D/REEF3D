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

#include"mooring_Catenary.h"
#include"lexer.h"
#include"ghostcell.h"
#include <Eigen/Dense>

mooring_Catenary::mooring_Catenary(int number):line(number){}

mooring_Catenary::~mooring_Catenary(){}

void mooring_Catenary::start(lexer *p, ghostcell *pgc)
{
    curr_time = p->simtime;

    FH_0 = 0.01;
    FV_0 = 0.01;

    calcForce(p,pgc);
		
	// Print mooring line
	print(p);
}


void mooring_Catenary::calcForce(lexer *p, ghostcell *pgc)
{
    Eigen::MatrixXd A_eigen = Eigen::MatrixXd::Zero(2,2); 
    Eigen::VectorXd B_eigen = Eigen::VectorXd::Zero(2);    
    Eigen::VectorXd F_eigen = Eigen::VectorXd::Zero(2);   

	// Calculate distances between start and mooring points
	dx = p->X311_xe[line] - p->X311_xs[line];			
	dy = p->X311_ye[line] - p->X311_ys[line];				
	dz = p->X311_ze[line] - p->X311_zs[line];	

	dxy_aim = sqrt(dx*dx+dy*dy);			

    double dxy_ = dxy_aim;

    for (int loop = 0; loop < 1000; loop++)
    {
        double f1, f2, df1H, df1V, df2H, df2V;
        FH = FH_0;
        FV = FV_0;
        FH_0 = FH + 1.0;
        FV_0 = FV + 1.0;	
       
        while (fabs(FH - FH_0) > 1e-5 && fabs(FV - FV_0) > 1e-5)
        {
            FH_0 = FH;
            FV_0 = FV;
   
            f1 = L - FV/w + FH/EA*FV/w + FH/w*log(FV/FH + sqrt(1+(FV/FH)*(FV/FH))) - dxy_;
            f2 = FH/w*(sqrt(1+(FV/FH)*(FV/FH))-1) + FV*FV/(2*EA*w) - dz;

            df1H = FV/(EA*w) + (FH*(-(FV/(FH*FH)) - (FV*FV)/((FH*FH*FH)*sqrt(1 + FV*FV/(FH*FH)))))/((FV/FH + sqrt(1 + FV*FV/(FH*FH)))*w) + log(FV/FH + sqrt(1 + FV*FV/(FH*FH)))/w;            
            df1V = (-1.0 + FH/EA)/w + (FH*((1/FH) + FV/(FH*FH*sqrt(1 + FV*FV/(FH*FH)))))/((FV/FH + sqrt(1 + FV*FV/(FH*FH)))*w);           
            df2H = (-1.0 + FH/sqrt(FH*FH + FV*FV))/w;
            df2V = FV/(EA*w) + FV/(FH*sqrt(1 + FV*FV/(FH*FH))*w);
           
            A_eigen << df1H, df1V, df2H, df2V;
            B_eigen << -f1, -f2;
            F_eigen = A_eigen.colPivHouseholderQr().solve(B_eigen);

            FH = FH + F_eigen[0];
            FV = FV + F_eigen[1];
        }
        Xme_ = FH*fabs(cos(atan(dy/dx)));
        Yme_ = FH*fabs(sin(atan(dy/dx)));
        Zme_ = FV;	

        buildLine(p);
   
        // Check convergence
        double dx_curr = x[H-1] - p->X311_xs[line];			
        double dy_curr = y[H-1] - p->X311_ys[line];			
        double dxy_curr = sqrt(dx_curr*dx_curr+dy_curr*dy_curr);			

        if (fabs(dxy_aim - dxy_curr) > 0.0001)
        {
            if (dxy_curr - dxy_aim > 0.0001)
            {
                dxy_ -= 0.0001;
            }
            else
            {
                dxy_ += 0.0001;
            }
        }
        else
        {
            break;
        }
    }

	// Reaction forces at mooring points	
	if (dx > 0)	
	{
		Xme_ *= -1.0;
	}
	if (dy > 0)	
	{
		Yme_ *= -1.0;
	}
	if (dz > 0)	
	{
		Zme_ *= -1.0;
	}
}


void mooring_Catenary::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
    // Tension forces if line is not broken
    if (broken == false)
    {
        Xme = Xme_; 
        Yme = Yme_;
        Zme = Zme_;
    }

    // Breakage due to max tension force
    if (breakTension > 0.0 && fabs(T[H-1]) >= breakTension)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }

    // Breakage due to time limit
    if (breakTime > 0.0 && curr_time >= breakTime)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }
}


void mooring_Catenary::getForce(lexer *p, ghostcell *pgc, double& FH_, double& FV_)
{
    // Ini line
	double rho_f = 1000.0;
	
	rho_c = p->X311_rho_c[line];
	w = p->X311_w[line]*9.81*(rho_c - rho_f)/rho_c;
	L = p->X311_l[line];
	H = p->X311_H[line];
	EA = p->X311_EA[line];
	
	p->Darray(x,H); 
	p->Darray(y,H);
	p->Darray(z,H); 
	p->Darray(T,H);

	printtime = 0.0;

    // Calculate shape
    FH_0 = 0.01;
    FV_0 = 0.01;
	calcForce(p, pgc);

    // Return values
	FH_ = FH;
	FV_ = FV;
}


void mooring_Catenary::getShape(lexer *p, ghostcell *pgc, double*& x_, double*& y_, double*& z_, double*& T_)
{
    // Ini line
	double rho_f = 1000.0;
	
	rho_c = p->X311_rho_c[line];
	w = p->X311_w[line]*9.81*(rho_c - rho_f)/rho_c;
	L = p->X311_l[line];
	H = p->X311_H[line] + 2;
	EA = p->X311_EA[line];
	
	p->Darray(x,H); 
	p->Darray(y,H);
	p->Darray(z,H); 
	p->Darray(T,H);

    // Calculate force
         
    double dx_ = x_[H-1] - x_[H-2];
    double dy_ = y_[H-1] - y_[H-2];
    double dxy_ = sqrt(dx_*dx_ + dy_*dy_);
    double dz_ = z_[H-1] - z_[H-2];

    double mag = sqrt(dxy_*dxy_ + dz_*dz_);

    FH_0 = T_[H-1]*dxy_/mag;
    FV_0 = T_[H-1]*dz_/mag;

    calcForce(p, pgc);

    // Return values
    x_ = x;
    y_ = y;
    z_ = z;
    T_ = T;
}


void mooring_Catenary::iniShape(lexer *p, ghostcell *pgc,Eigen::VectorXd& x_, Eigen::VectorXd& y_, Eigen::VectorXd& z_)
{
    // Ini line
	double rho_f = 1000.0;
	
	rho_c = p->X311_rho_c[line];
	w = p->X311_w[line]*9.81*(rho_c - rho_f)/rho_c;
	L = p->X311_l[line];
	H = p->X311_H[line] + 1;
	EA = p->X311_EA[line];
	
    p->Darray(x,H); 
	p->Darray(y,H);
	p->Darray(z,H); 
	p->Darray(T,H);

    // Calculate force
         
    FH_0 = 1.0;
    FV_0 = 1.0; 

    calcForce(p, pgc);

    // Return values
    for (int ii = 0; ii < H; ii++)
    {
        x_(ii) = x[ii];
        y_(ii) = y[ii];
        z_(ii) = z[ii];
    }
}
