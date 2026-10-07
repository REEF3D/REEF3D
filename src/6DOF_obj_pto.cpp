/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

// PTO of the FNPF body (X 500 - X 504): joint, element set-up, loads, output.
//
// X 500 1: explicit, the PTO load enters F like the mooring load
// X 500 2: the Jacobians of the stiff elements (end-stops) are added to the added-mass system
//          (M + A - h dF_s/dqd J J^T - h^2 dF_s/dq J J^T) a = F + h dF_s/dq qd J,  h = dt
//          linearly implicit, stable for stiff end-stops, first order during contact;
//          dampers and springs stay explicit (RK3/RK4 accuracy)

#include"6DOF_obj.h"
#include"6DOF_pto_elements.h"
#include"6DOF_output_dir.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>

void sixdof_obj::ini_pto(lexer *p, ghostcell *pgc)
{
    pto_on_ = (p->X500>0);
    pto_implicit_ = (p->X500==2);

    if(!pto_on_)
    return;

    // joint: axis and attachment point in world coordinates, default heave at the CoG
    Eigen::Vector3d axis(p->X501_ax, p->X501_ay, p->X501_az);
    Eigen::Vector3d xatt = c_;

    if(axis.norm()<1.0e-12)
    {
        if(p->mpirank==0)
        cout<<"6DOF PTO: X 501 axis has zero length"<<endl;

        pgc->final();
        exit(1);
    }

    if(p->X501==1)
    {
        // world -> model, the axis as the difference of two transformed points
        double x1=p->X501_xa, y1=p->X501_ya;
        double x2=p->X501_xa+axis(0), y2=p->X501_ya+axis(1);

        p->XYin(x1,y1);
        p->XYin(x2,y2);

        xatt = Eigen::Vector3d(x1, y1, p->X501_za);
        axis = Eigen::Vector3d(x2-x1, y2-y1, axis(2));
    }

    pto_joint_.initialize(axis,xatt,c_,R_);

    // elements
    if(p->X502==1)
    pto_.add(std::make_unique<pto_damper>(p->X502_B));

    if(p->X503==1)
    pto_.add(std::make_unique<pto_springdamper>(p->X503_K, p->X503_B));

    if(p->X504==1)
    pto_.add(std::make_unique<pto_endstop>(p->X504_qmin, p->X504_qmax, p->X504_K, p->X504_C));

    if(pto_.empty() && p->mpirank==0)
    cout<<"6DOF PTO: X 500 set but no PTO element (X 502 - X 504) given"<<endl;

    pto_E_ = 0.0;
    pto_Pn_ = 0.0;
    pto_tn_ = -1.0;
    pto_warned_ = false;

    if(p->mpirank==0)
    {
        cout<<"6DOF PTO: axis "<<pto_joint_.a.transpose()<<"  attachment "<<xatt.transpose()
        <<"  "<<(pto_implicit_?"implicit":"explicit")<<endl;

        for(auto &e : pto_.elem)
        cout<<"6DOF PTO element: "<<e->name()<<endl;

        char str[1000];
        sprintf(str,"%s/REEF3D_6DOF_pto_%i.dat",sixdof_output_dir(p),n6DOF);

        printpto.open(str);
        printpto<<"time \t q \t qd \t F_pto \t P_abs \t E_abs";

        for(auto &e : pto_.elem)
        printpto<<" \t F_"<<e->name();

        printpto<<endl;
    }
}

void sixdof_obj::pto_forces(lexer *p, ghostcell*, int iter)
{
    // state of the current RK stage (u_fb: also prescribed motions, X 11 2)
    const Eigen::Vector3d v(u_fb(0), u_fb(1), u_fb(2));
    const Eigen::Vector3d w(u_fb(3), u_fb(4), u_fb(5));

    pto_s_ = pto_joint_.state(c_,R_,v,w,p->simtime);
    pto_out_ = pto_.force(pto_s_);

    const Eigen::Matrix<double,6,1> Fg = pto_out_.F*pto_joint_.J;

    Xext += Fg(0);
    Yext += Fg(1);
    Zext += Fg(2);
    Kext += Fg(3);
    Mext += Fg(4);
    Next += Fg(5);

    // first stage: accepted state at t^n
    if(iter==0)
    {
        if(pto_tn_>=0.0 && p->simtime>pto_tn_)
        pto_E_ += 0.5*(pto_out_.P + pto_Pn_)*(p->simtime - pto_tn_);

        pto_Pn_ = pto_out_.P;
        pto_tn_ = p->simtime;

        pto_.accept(pto_s_,p->dt);

        pto_stability_check(p);

        if(p->mpirank==0)
        {
            printpto<<p->simtime<<" \t "<<pto_s_.q<<" \t "<<pto_s_.qd<<" \t "<<pto_out_.F
            <<" \t "<<pto_out_.P<<" \t "<<pto_E_;

            for(auto &e : pto_.elem)
            printpto<<" \t "<<e->F_last;

            printpto<<endl;

            if(p->count%p->P12==0)
            cout<<"PTO q: "<<pto_s_.q<<" qd: "<<pto_s_.qd<<" F: "<<pto_out_.F<<" P: "<<pto_out_.P<<" E: "<<pto_E_<<endl;
        }
    }
}

void sixdof_obj::pto_implicit(lexer *p, Eigen::Matrix<double,6,6> &L, Eigen::Matrix<double,6,1> &r)
{
    if(!(pto_on_ && pto_implicit_))
    return;

    const double h = p->dt;
    const Eigen::Matrix<double,6,1> &J = pto_joint_.J;

    L -= (h*pto_out_.dFdqd_s + h*h*pto_out_.dFdq_s)*(J*J.transpose());
    r += (h*pto_out_.dFdq_s*pto_s_.qd)*J;
}

void sixdof_obj::pto_stability_check(lexer *p)
{
    // explicit RK3/RK4: |lambda dt| below ~2.5 for the explicitly treated PTO damping and
    // stiffness, with the effective mass at the joint m_eff = 1/(J^T (M + A)^-1 J) of the free DOFs
    if(pto_warned_ || !am_on_)
    return;

    Eigen::Matrix<double,6,6> L = Aadd_;
    L.block<3,3>(0,0) += Mass_fb*Eigen::Matrix3d::Identity();
    L.block<3,3>(3,3) += R_*I_*R_.transpose();

    Eigen::Matrix<double,6,1> J = pto_joint_.J;

    for(int n=0; n<6; ++n)
    if(p_fixed_dof(p,n))
    {
    L.row(n).setZero();
    L.col(n).setZero();
    L(n,n) = 1.0;
    J(n) = 0.0;
    }

    const double JLJ = J.dot(L.partialPivLu().solve(J));

    if(JLJ<=1.0e-20)
    return;

    const double meff = 1.0/JLJ;
    const double Kx = pto_out_.dFdq  - (pto_implicit_ ? pto_out_.dFdq_s  : 0.0);
    const double Bx = pto_out_.dFdqd - (pto_implicit_ ? pto_out_.dFdqd_s : 0.0);

    const double lam = fabs(Bx)/meff;
    const double om = sqrt(fabs(Kx)/meff);
    const double cfl = p->dt*MAX(lam,om);

    if(cfl>2.5)
    {
        pto_warned_ = true;

        if(p->mpirank==0)
        {
        cout<<"6DOF PTO: explicit coupling at its stability limit (dt*lambda = "<<cfl<<")";
        cout<<(pto_implicit_ ? ", reduce the time step" : ", use X 500 2 or reduce the time step")<<endl;
        }
    }
}
