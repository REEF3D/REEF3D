/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2021 Tobias Martin

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

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsi_strip::setFieldBC(Matrix3Xd& c_, Matrix3Xd& cdot_, Matrix4Xd& q_, Matrix4Xd& q0_, Matrix4Xd& qdot_, Matrix3Xd& f_, Matrix4Xd& m0_, Matrix3Xd& rhs_cdot_, double time , int ind)
{
    if (ind == 0)
    {
        //BC: Fixed rotatory end
        Eigen::Vector4d qb; 
        Eigen::Matrix4d J; 

        qb << 1.0, 0, 0, 0;
        J << 2.0*qb(0)*qb(0) - 1.0, 2*qb(1)*qb(0), 2*qb(2)*qb(0), 2*qb(3)*qb(0),
             2*qb(0)*qb(1), 2.0*qb(1)*qb(1) - 1.0, 2*qb(2)*qb(1), 2*qb(3)*qb(1),
             2*qb(0)*qb(2), 2*qb(1)*qb(2), 2.0*qb(2)*qb(2) - 1.0, 2*qb(3)*qb(2),
             2*qb(0)*qb(3), 2*qb(1)*qb(3), 2*qb(2)*qb(3), 2.0*qb(3)*qb(3) - 1.0;
         
        q_.col(0) = 2.0*qb.dot(q_.col(1))*qb - q_.col(1);
        qdot_.col(0) = J*qdot_.col(1); 
    }
    else if (ind == 1)
    {
        // BC: Free translatory end with vanishing forces
        f_.col(Ne+1) = -f_.col(Ne); // correct?
    }
    else if (ind == 2)
    {
        // BC: Free rotatory end with vanishing moments 
        m0_.col(Ne) = Eigen::Vector4d::Zero(4); q_.col(Ne+1) = q_.col(Ne); q0_.col(Ne+1) = q0_.col(Ne);
    }
    else if (ind == 3)
    {
        // BC: Fixed translatory end
        rhs_cdot_.col(0) = Eigen::Vector3d::Zero(3);
    }
}


void fsi_strip::setConstantLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_)
{
}


void fsi_strip::setVariableLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_, const double time)
{
    double dm_el, m_el;
    Eigen::Vector3d P_el_star, I_el_star, s0, omega_el, omega_el_0;
    Eigen::Matrix3d J0, Xil_0_skew;

    double delta_time = (time - t_strip_n)/(t_strip - t_strip_n) > 0.01 ? time - t_strip_n : 1e20;

    for (int eI = 1; eI < Ne+1; eI++)
    {
        m_el = 0.0;   
        P_el.col(eI) << 0.0, 0.0, 0.0;
        P_el_star << 0.0, 0.0, 0.0;
        I_el_star << 0.0, 0.0, 0.0;
        s0 << 0.0, 0.0, 0.0;
        J0 << Eigen::Matrix3d::Zero();
        
        for (int pI = 0; pI < lagrangePoints[eI-1].cols(); pI++)
        {
            // Mass of element
            dm_el = rho_f*dx_body*lagrangeArea[eI-1](pI);
            m_el += dm_el;

            // Preliminary linear momentum
            P_el_star += dm_el*lagrangeVel[eI-1].col(pI);

            // Static moment
            s0 += dm_el*Xil_0[eI-1].col(pI);
            
            // Preliminary angular momentum
            I_el_star += dm_el*Xil[eI-1].col(pI).cross(lagrangeVel[eI-1].col(pI));

            // Quaternionic tensor of inertia
            Xil_0_skew << 0, -Xil_0[eI-1](2,pI), Xil_0[eI-1](1,pI), Xil_0[eI-1](2,pI), 0, -Xil_0[eI-1](0,pI), -Xil_0[eI-1](1,pI), Xil_0[eI-1](0,pI), 0;
            Xil_0_skew = Xil_0_skew.transpose()*Xil_0_skew;
            J0 += dm_el*Xil_0_skew;
        }

        // Determine linear momentum
        omega_el = getOmega(q_.col(eI),qdot_.col(eI));
        P_el.col(eI) = m_el*(cdot_.col(eI-1)+cdot_.col(eI))/2.0 + omega_el.cross(rotVec(s0,q_.col(eI)));

        // Determine coupling force
        F_el.col(eI) = -(P_el.col(eI) - P_el_n.col(eI))/delta_time - (P_el_n.col(eI) - P_el_star)/(t_strip - t_strip_n);

        // Determine angular momentum
        omega_el_0 = getOmega0(q_.col(eI),qdot_.col(eI));
        I_el.col(eI) = rotVec(s0,q_.col(eI)).cross((cdot_.col(eI-1)+cdot_.col(eI))/2.0) + rotVec(J0*omega_el_0,q_.col(eI));
        
        // Determine coupling moment
        M_el.col(eI) = -(I_el.col(eI) - I_el_n.col(eI))/delta_time - (I_el_n.col(eI) - I_el_star)/(t_strip - t_strip_n);
    }
    
    // Assign external forces
    for (int eI = 0; eI < Ne+1; eI++)
    {
        Fext_.col(eI) = (1.0 - rho_f/rho_s)*gravity_vec + (F_el.col(eI) + F_el.col(eI+1))/(2.0*rho_s*A_el*l_el);
    }
    
    // Assign external moments
    for (int eI = 0; eI < Ne+2; eI++)
    {
        Mext_.col(eI) << 0.0, M_el.col(eI)/l_el;
    }
}
