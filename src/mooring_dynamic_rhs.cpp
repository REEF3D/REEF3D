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

#include"mooring_dynamic.h"
#include"lexer.h"
#include"ghostcell.h"

void mooring_dynamic::setConstantLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_)
{
    double zg = 0.0;
    double cd_t = 0.5;
    double cd_n = 2.5;
    double cm_t = 0.0;
    double cm_n = 3.8;
    double Kg = 3.0e6;
    double nu = 0.01;
    double mu = 0.3;
    double xi = 1.0;

    double t_x, t_y, t_z, t_mag, v_x, v_y, v_z, vn_x, vn_y, vn_z, vn_mag, v_t;
    double a_t, a_x, a_y, a_z, an_x, an_y, an_z;
        
    double Fb = -9.81*(rho_c - 1000)/rho_c;

    for (int i = 0; i < Ne+1; i++)
    {
        // Direction vectors
        if (i < Ne)
        {
            t_x = c_(0,i+1) - c_(0,i); 
            t_y = c_(1,i+1) - c_(1,i); 
            t_z = c_(2,i+1) - c_(2,i); 
        }
        else
        {
            t_x = c_(0,i) - c_(0,i-1); 
            t_y = c_(1,i) - c_(1,i-1); 
            t_z = c_(2,i) - c_(2,i-1); 
        }
        t_mag = sqrt(t_x*t_x + t_y*t_y + t_z*t_z);
        t_x = t_x/t_mag; t_y = t_y/t_mag; t_z = t_z/t_mag;
        
        // Velocity vectors
        v_x = fluid_vel[i][0] - cdot_(0,i);
        v_y = fluid_vel[i][1] - cdot_(1,i);
        v_z = fluid_vel[i][2] - cdot_(2,i);

        // Acceleration vectors
        a_x = fluid_acc[i][0] - cdotdot_moor(0,i);
        a_y = fluid_acc[i][1] - cdotdot_moor(1,i);
        a_z = fluid_acc[i][2] - cdotdot_moor(2,i);
        
        // Drag force
        v_t = v_x*t_x + v_y*t_y + v_z*t_z;
        vn_x = v_x - v_t*t_x;
        vn_y = v_y - v_t*t_y;
        vn_z = v_z - v_t*t_z;
        vn_mag = sqrt(vn_x*vn_x + vn_y*vn_y + vn_z*vn_z);
    
        Fext_(0,i) = 0.5*rho_c*d_c*(cd_t*fabs(v_t)*v_t*t_x + cd_n*vn_mag*vn_x);
        Fext_(1,i) = 0.5*rho_c*d_c*(cd_t*fabs(v_t)*v_t*t_y + cd_n*vn_mag*vn_y);
        Fext_(2,i) = 0.5*rho_c*d_c*(cd_t*fabs(v_t)*v_t*t_z + cd_n*vn_mag*vn_z);

        // Added mass force
        a_t = a_x*t_x+a_y*t_y+a_z*t_z;
        an_x = a_x - a_t*t_x;
        an_y = a_y - a_t*t_y;
        an_z = a_z - a_t*t_z;
    
        Fext_(0,i) += rho_c*PI/4.0*d_c*d_c*(fluid_acc[i][0] + cm_t*a_t*t_x + cm_n*an_x);
        Fext_(1,i) += rho_c*PI/4.0*d_c*d_c*(fluid_acc[i][1] + cm_t*a_t*t_y + cm_n*an_y);
        Fext_(2,i) += rho_c*PI/4.0*d_c*d_c*(fluid_acc[i][2] + cm_t*a_t*t_z + cm_n*an_z);

        // Gravity force
        Fext_(2,i) += Fb;

        // Bottom force
        if ((zg - c_(2,i)) > 0.0)
        {
            Fext_(2,i) += (Kg*d_c*(zg - c_(2,i)) - 2.0*xi*sqrt(Kg*rho_c*A*d_c)*max(cdot_(2,i),0.0));
            
            double vx = cdot_(0,i)/max(nu,sqrt(cdot_(0,i)*cdot_(0,i) + cdot_(1,i)*cdot_(1,i)));
            double vy = cdot_(1,i)/max(nu,sqrt(cdot_(0,i)*cdot_(0,i) + cdot_(1,i)*cdot_(1,i)));
            Fext_(0,i) = mu*Fb*sin(PI/2.0*vx);
            Fext_(1,i) = mu*Fb*sin(PI/2.0*vy);
        }
    }
}

void mooring_dynamic::setFieldBC(Matrix3Xd& c_, Matrix3Xd& cdot_, Matrix4Xd& q_, Matrix4Xd& q0_, Matrix4Xd& qdot_, Matrix3Xd& f_, Matrix4Xd& m0_, Matrix3Xd& rhs_cdot_, double time , int ind)
{
    if (ind==0)
    {
        double delta_tm = t_mooring - t_mooring_n;
        double tau = time - t_mooring_n;

        // Quadratic interpolation
        //a_O = 1.0/(1.0/2.0*delta_tm*delta_tm)*(fixPoint - c_moor_n.col(Ne) - cdot_n.col(Ne)*delta_tm);
        //cdot_.col(Ne) = cdot_n.col(Ne) + a_O*tau;
        //c_.col(Ne) = c_moor_n.col(Ne) + cdot_n.col(Ne)*tau + 0.5*a_O*tau*tau;

        // Linear interpolation in solver
        cdot_.col(Ne) = (fixPoint - c_moor_n.col(Ne))/delta_tm;
    }
    else if (ind==1)
    {
        // BC: Free translatory end with vanishing forces
        // f.col(Ne+1) = -f.col(Ne); // correct?
    }
    else if (ind==2)
    {
        // BC: Free rotatory end with vanishing moments 
        m0_.col(0) = Eigen::Vector4d::Zero(4); q_.col(0) = q_.col(1); q0_.col(0) = q0_.col(1); 
        m0_.col(Ne) = Eigen::Vector4d::Zero(4); q_.col(Ne+1) = q_.col(Ne); q0_.col(Ne+1) = q0_.col(Ne);
    }
    else if (ind==3)
    {
        // BC: Fixed translatory end
        rhs_cdot_.col(0) = Eigen::Vector3d::Zero(3);
        rhs_cdot_.col(Ne) = Eigen::Vector3d::Zero(3);
    }
}


