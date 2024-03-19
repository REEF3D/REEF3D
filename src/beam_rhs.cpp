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

#include"beam.h"

void beam::rhs(Matrix3Xd& c_, Matrix3Xd& cdot_, Matrix4Xd& q_, Matrix4Xd& qdot_, double time)
{
    // Set variable loads 
    setVariableLoads(Fext, Mext, c_, cdot_, q_, qdot_, time);

    // Set field boundary conditions
    setFieldBC(c_, cdot_, q_, q0, qdot_, f, m0, rhs_cdot, time, 0);
    
    // Internal forces
    for (int i = 1; i < Ne+1; i++)
    {
        f0.col(i) = f0_(c_.col(i),c_.col(i-1),c0.col(i),c0.col(i-1),cdot_.col(i),cdot_.col(i-1),q_.col(i),qdot_.col(i),q0.col(i));
        f.col(i) = R(q_.col(i))*f0.col(i).tail(3);
    }

    // Set field boundary conditions
    setFieldBC(c_, cdot_, q_, q0, qdot_, f, m0, rhs_cdot, time, 1);
    
    // Internal moments
    for (int i = 0; i < Ne+1; i++)
    {   
        m0.col(i) = m0_(q_.col(i+1),q_.col(i),q0.col(i+1),q0.col(i),qdot_.col(i+1),qdot_.col(i));
    }
    
    // Set field boundary conditions
    setFieldBC(c_, cdot_, q_, q0, qdot_, f, m0, rhs_cdot, time, 2);

    // Determine internal RHS
    rhs_cdot = 
        1.0/(rho*A)*
        (
            (f.block(0,1,3,Ne+1) - f.block(0,0,3,Ne+1))/dZ
        )
        + Fext; 
   
    for (int i = 1; i < Ne+1; i++)
    {
        calcInvM(q_.col(i));    

        rhs_qdot.col(i) = 
            2.0/rho*invM*
            ( 
                4.0*rho*qMult(qdot_.col(i), Iq*qMult(qconj(qdot_.col(i)),q_.col(i)))
              + qMult(0, (c_.col(i) - c_.col(i-1))/dZ, q_.col(i), f0.col(i)) 
              + 1.0/dZ*(qMult(q_.col(i+1),m0.col(i)) - qMult(q_.col(i-1),m0.col(i-1))) 
              + qMult(Mext.col(i),q_.col(i))
            )
            - qdot_.col(i).dot(qdot_.col(i))*q_.col(i);
    }

    // Set field boundary conditions
    setFieldBC(c_, cdot_, q_, q0, qdot_, f, m0, rhs_cdot, time, 3);
}


Eigen::Vector4d beam::f0_
(
    const Eigen::Vector3d& cr, const Eigen::Vector3d& cl,
    const Eigen::Vector3d& c0r, const Eigen::Vector3d& c0l,
    const Eigen::Vector3d& cdotr, const Eigen::Vector3d& cdotl,
    const Eigen::Vector4d& qi, const Eigen::Vector4d& qdoti,
    const Eigen::Vector4d& q0i
)
{
    dcdz.tail(3)  = (cr - cl)/dZ;  
    dc0dz.tail(3) = (c0r - c0l)/dZ;
    dcdotdz.tail(3) = (cdotr - cdotl)/dZ;

    fdot = qMult(qconj(qdoti),dcdz,qi) + qMult(qconj(qi),dcdotdz,qi) + qMult(qconj(qi),dcdz,qdoti);
  
    Eigen::Vector4d force = Eigen::Vector4d::Zero(4);
    
    force.tail(3) = Ceps*(R(qi).transpose()*dcdz.tail(3) - R(q0i).transpose()*dc0dz.tail(3)) + 2.0*Cepsdot*fdot.tail(3);

    // f0 calculation without compression effects 
    if (compression==false)
    {
        Eigen::Vector3d f0_ini; f0_ini << 1,0,0;
        Eigen::Vector3d f0_cur = R(qi).transpose()*dcdz.tail(3);
        if (f0_cur(0)-f0_ini(0) < 0.0)
        {
            f0_cur(0) = 0.0;
            f0_ini(0) = 0.0;
        }
        force.tail(3) = Ceps*(f0_cur - f0_ini) + 2.0*Cepsdot*fdot.tail(3);
    }

    return force;
}


    Eigen::Vector4d beam::m0_
(
    const Eigen::Vector4d& qr, const Eigen::Vector4d& ql,
    const Eigen::Vector4d& q0r, const Eigen::Vector4d& q0l,
    const Eigen::Vector4d& qdotr, const Eigen::Vector4d& qdotl
)
{
    Eigen::Vector4d moment = Eigen::Vector4d::Zero(4);

    // Strain
    Eigen::Vector4d mult = qMult(qconj(ql),qr);
    Eigen::Vector3d kappa = 2.0/dZ*sqrt(2.0/(1.0 + mult(0)))*mult.tail(3);
    
    moment.tail(3) = Ckappa*kappa;
    
    // Initial strain
    mult = qMult(qconj(q0l),q0r);
    kappa = 2.0/dZ*sqrt(2.0/(1.0 + mult(0)))*mult.tail(3);

    moment.tail(3) -= Ckappa*kappa;

    // Strain velocity
    /*mult = qMult(qconj(qdotl),qr);
    kappa = 2.0/dZ*sqrt(2.0/(1.0 + mult(0)))*mult.tail(3);

    mult = qMult(qconj(ql),qdotr);
    kappa += 2.0/dZ*sqrt(2.0/(1.0 + mult(0)))*mult.tail(3);
    */

    mult = qMult(qconj(ql),qr);
    double zeta = sqrt(2.0/(1.0 + mult(0)));
    dummy = qMult(qconj(qdotl),qr) + qMult(qconj(ql),qdotr);
    double zeta_ = -0.25*pow(zeta,3);            
    kappa = 2.0/dZ*(zeta_*dummy(0)*mult.tail(3) + zeta*dummy.tail(3));

    moment.tail(3) += 2.0*Ckappadot*kappa;

    return moment;
}

void beam::calcInvM(const Eigen::Vector4d& q)
{
    calcQ(q);
    invM.noalias() = 0.25*Q*invIq*Q.transpose();     
}

Eigen::Vector4d beam::qconj(const Eigen::Vector4d& q)
{
    dummy << q(0), -q(1), -q(2), -q(3);
    return dummy;
}

void beam::calcQ(const Eigen::Vector4d& q)
{
    Q << q(0), -q(1), -q(2), -q(3),
         q(1),  q(0), -q(3),  q(2),
         q(2),  q(3),  q(0), -q(1),
         q(3), -q(2),  q(1),  q(0);
}
    
Eigen::Vector4d beam::qMult(const Eigen::Vector4d& q1, const Eigen::Vector4d& q2)
{
    calcQ(q1);
    dummy.noalias() = Q*q2;
    return dummy;
}
    
Eigen::Vector4d beam::qMult(const Eigen::Vector4d& q1, const Eigen::Vector4d& q2, const Eigen::Vector4d& q3)
{
    Eigen::Vector4d q12;
    calcQ(q1);
    q12.noalias() = Q*q2;
    calcQ(q12);
    dummy.noalias() = Q*q3;
    return dummy;
}

Eigen::Vector4d beam::qMult(int, const Eigen::Vector3d& c1, const Eigen::Vector4d& q2, const Eigen::Vector4d& q3)
{
    dummy = Eigen::Vector4d::Zero(4);
    dummy.tail(3) = c1;
    calcQ(dummy);
    Eigen::Vector4d q12;
    q12.noalias() = Q*q2;
    calcQ(q12);
    dummy.noalias() = Q*q3;
    return dummy;
}

Eigen::Matrix3d beam::R(const Eigen::Vector4d& e)
{
    Eigen::Matrix<double,3,4> E,G;
    Eigen::Matrix3d R_;

    E << -e(1), e(0), -e(3), e(2),
         -e(2), e(3), e(0), -e(1),
         -e(3), -e(2), e(1), e(0); 

    G << -e(1), e(0), e(3), -e(2),
         -e(2), -e(3), e(0), e(1),
         -e(3), e(2), -e(1), e(0); 

    R_.noalias() = E*G.transpose(); 

    return R_; 
}

Eigen::Vector3d beam::getOmega(const Eigen::Vector4d& qI, const Eigen::Vector4d& qdotI)
{
    return 2.0*qMult(qdotI,qconj(qI)).tail(3);
}

Eigen::Vector3d beam::getOmega0(const Eigen::Vector4d& qI, const Eigen::Vector4d& qdotI)
{
    return 2.0*qMult(qconj(qI),qdotI).tail(3);
}
    
Eigen::Vector3d beam::rotVec(const Eigen::Vector3d& vec_, const Eigen::Vector4d& qI)
{
    return R(qI)*vec_;
}

