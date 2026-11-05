/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsi_strip::setFieldBC(Matrix3Xd& c_, Matrix3Xd& cdot_, Matrix4Xd& q_, Matrix4Xd& q0_, Matrix4Xd& qdot_, Matrix3Xd& f_, Matrix4Xd& m0_, Matrix3Xd& rhs_cdot_, double time , int ind)
{
    if (ind==0)
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
    else if (ind==1)
    {
        // BC: Free translatory end with vanishing forces
        f_.col(Ne+1) = -f_.col(Ne); // correct?
    }
    else if (ind==2)
    {
        // BC: Free rotatory end with vanishing moments 
        m0_.col(Ne) = Eigen::Vector4d::Zero(4); q_.col(Ne+1) = q_.col(Ne); q0_.col(Ne+1) = q0_.col(Ne);
    }
    else if (ind==3)
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
    // Called in every RHS evaluation of the beam solver. The element mass, static moment and
    // inertia (m_el, s0, J0) are constant, and the fluid momentum (P_el_star, I_el_star) does
    // not change during one Integrate call, so the loops over the Lagrangian points were moved
    // to precompute_element_constants() and precompute_fluid_momentum().
    
    double m_el;
    Eigen::Vector3d P_el_star, I_el_star, s0, omega_el, omega_el_0;
    Eigen::Matrix3d J0;

    double delta_time = (time - t_strip_n)/(t_strip - t_strip_n) > 0.01 ? time - t_strip_n : 1e20;

    for (int eI = 1; eI < Ne+1; eI++)
    {
        m_el = m_el_c(eI-1);
        s0 = s0_c.col(eI-1);
        J0 = J0_c[eI-1];
        P_el_star = P_star_c.col(eI-1);
        I_el_star = I_star_c.col(eI-1);

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
    if (thinStrip==false)
    {
        for (int eI = 0; eI < Ne+2; eI++)
        {
            Mext_.col(eI) << 0.0, M_el.col(eI)/l_el;
        }
    }
}

void fsi_strip::precompute_element_constants()
{
    double dm_el;
    Eigen::Matrix3d Xil_0_skew;
    
    m_el_c = Eigen::VectorXd::Zero(Ne);
    s0_c = Matrix3Xd::Zero(3,Ne);
    J0_c.assign(Ne, Eigen::Matrix3d::Zero());
    
    for (int eI = 0; eI < Ne; eI++)
    {
        double m_el = 0.0;
        Eigen::Vector3d s0 = Eigen::Vector3d::Zero();
        Eigen::Matrix3d J0 = Eigen::Matrix3d::Zero();
        
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            // Mass of element
            dm_el = rho_f*dx_body*lagrangeArea[eI](pI);
            m_el += dm_el;

            // Static moment
            s0 += dm_el*Xil_0[eI].col(pI);

            // Quaternionic tensor of inertia
            Xil_0_skew << 0, -Xil_0[eI](2,pI), Xil_0[eI](1,pI), Xil_0[eI](2,pI), 0, -Xil_0[eI](0,pI), -Xil_0[eI](1,pI), Xil_0[eI](0,pI), 0;
            Xil_0_skew = Xil_0_skew.transpose()*Xil_0_skew;
            J0 += dm_el*Xil_0_skew;
        }
        
        m_el_c(eI) = m_el;
        s0_c.col(eI) = s0;
        J0_c[eI] = J0;
    }
    
    P_star_c = Matrix3Xd::Zero(3,Ne);
    I_star_c = Matrix3Xd::Zero(3,Ne);
}

void fsi_strip::precompute_fluid_momentum()
{
    double dm_el;
    
    for (int eI = 0; eI < Ne; eI++)
    {
        Eigen::Vector3d P_el_star = Eigen::Vector3d::Zero();
        Eigen::Vector3d I_el_star = Eigen::Vector3d::Zero();
        
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            dm_el = rho_f*dx_body*lagrangeArea[eI](pI);

            // Preliminary linear momentum
            P_el_star += dm_el*lagrangeVel[eI].col(pI);
            
            // Preliminary angular momentum
            I_el_star += dm_el*Xil[eI].col(pI).cross(lagrangeVel[eI].col(pI));
        }
        
        P_star_c.col(eI) = P_el_star;
        I_star_c.col(eI) = I_el_star;
    }
}
