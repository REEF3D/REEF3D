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

// PTO of the FNPF body (X 500 - X 508): joint, element set-up, controllers, loads, output.
//
// X 500 1: explicit, the PTO load enters F like the mooring load
// X 500 2: the Jacobians of the stiff elements (end-stops) are added to the added-mass system
//          (M + A - h dF_s/dqd J J^T - h^2 dF_s/dq J J^T) a = F + h dF_s/dq qd J,  h = dt
//          linearly implicit, stable for stiff end-stops, first order during contact;
//          dampers and springs stay explicit (RK3/RK4 accuracy)
// X 505: generator force and power limits, X 506: tuned passive damping (pto_hydro.dat),
// X 507: latching, X 508: declutching

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

    if(p->X507==1)
    pto_.add(std::make_unique<pto_latch>(p->X507_K, p->X507_C));

    if(p->X504==1)
    pto_.add(std::make_unique<pto_endstop>(p->X504_qmin, p->X504_qmax, p->X504_K, p->X504_C));

    // generator limits
    if(p->X505==1)
    {
    pto_.Fmax = p->X505_Fmax;
    pto_.Pmax = p->X505_Pmax;
    }

    // controllers
    if(p->X506==1)
    pto_tuned_passive(p,pgc);

    if(p->X507==1)
    pto_ctrl_.push_back(std::make_unique<pto_ctrl_latching>(p->X507_t));

    if(p->X508==1)
    pto_ctrl_.push_back(std::make_unique<pto_ctrl_declutching>(p->X508_t));

    if(pto_.empty() && p->mpirank==0)
    cout<<"6DOF PTO: X 500 set but no PTO element (X 502 - X 504, X 506, X 507) given"<<endl;

    if(p->X507==1 && !pto_implicit_ && p->mpirank==0)
    cout<<"6DOF PTO: latching (X 507) with explicit coupling, X 500 2 recommended for a stiff latch"<<endl;

    pto_E_ = 0.0;
    pto_Estep_ = 0.0;
    pto_warned_ = false;

    if(p->mpirank==0)
    {
        cout<<"6DOF PTO: axis "<<pto_joint_.a.transpose()<<"  attachment "<<xatt.transpose()
        <<"  "<<(pto_implicit_?"implicit":"explicit")<<endl;

        for(auto &e : pto_.elem)
        cout<<"6DOF PTO element: "<<e->name()<<endl;

        for(auto &c : pto_ctrl_)
        cout<<"6DOF PTO controller: "<<c->name()<<endl;

        if(pto_.Fmax>0.0 || pto_.Pmax>0.0)
        cout<<"6DOF PTO limits: Fmax "<<pto_.Fmax<<" N  Pmax "<<pto_.Pmax<<" W"<<endl;

        char str[1000];
        sprintf(str,"%s/REEF3D_6DOF_pto_%i.dat",sixdof_output_dir(p),n6DOF);

        printpto.open(str);
        printpto<<"time \t q \t qd \t F_pto \t P_abs \t E_abs \t gain \t latch \t sat";

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

    // stage time: TVD RK3 c = 0, 1, 1/2; RK4 c = 0, 1/2, 1/2, 1 (p->simtime is t^n in all stages)
    static const double crk3[3] = {0.0, 1.0, 0.5};
    static const double crk4[4] = {0.0, 0.5, 0.5, 1.0};
    double ts = p->simtime;

    if(p->A310==3 && iter>=0 && iter<3)
    ts += crk3[iter]*p->dt;

    else if(p->A310!=3 && iter>=0 && iter<4)
    ts += crk4[iter]*p->dt;

    pto_s_ = pto_joint_.state(c_,R_,v,w,ts);

    // controllers act once per time step on the accepted state at t^n
    if(iter==0)
    {
        pto_ctrl_in in;
        in.t  = p->simtime;
        in.q  = pto_s_.q;
        in.qd = pto_s_.qd;

        for(auto &ctrl : pto_ctrl_)
        ctrl->update(in,pto_);
    }

    pto_out_ = pto_.force(pto_s_);

    const Eigen::Matrix<double,6,1> Fg = pto_out_.F*pto_joint_.J;

    Xext += Fg(0);
    Yext += Fg(1);
    Zext += Fg(2);
    Kext += Fg(3);
    Mext += Fg(4);
    Next += Fg(5);

    // absorbed energy with the weights of the RK scheme of the body (TVD RK3: 1/6, 1/6, 2/3;
    // RK4: 1/6, 1/3, 1/3, 1/6), i.e. the work the integrator applies, also across switches
    if(iter==0)
    {
    pto_E_ += pto_Estep_;
    pto_Estep_ = 0.0;
    }

    static const double wrk3[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    static const double wrk4[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

    if(p->A310==3 && iter>=0 && iter<3)
    pto_Estep_ += wrk3[iter]*pto_out_.P*p->dt;

    else if(p->A310!=3 && iter>=0 && iter<4)
    pto_Estep_ += wrk4[iter]*pto_out_.P*p->dt;

    // first stage: accepted state at t^n
    if(iter==0)
    {
        pto_.accept(pto_s_,p->dt);

        pto_stability_check(p);

        if(p->mpirank==0)
        {
            double latch=0.0;

            for(auto &ctrl : pto_ctrl_)
            if(std::string(ctrl->name())=="latching")
            latch = ctrl->state();

            const double gain = (p->simtime>=pto_.t_off0 && p->simtime<pto_.t_off1) ? 0.0 : pto_.gain;

            printpto<<p->simtime<<" \t "<<pto_s_.q<<" \t "<<pto_s_.qd<<" \t "<<pto_out_.F
            <<" \t "<<pto_out_.P<<" \t "<<pto_E_<<" \t "<<gain<<" \t "<<latch<<" \t "<<pto_.sat_last;

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

void sixdof_obj::pto_tuned_passive(lexer *p, ghostcell *pgc)
{
    // B_opt at the period T (X 506; T <= 0: the wave input period B 93, Tp for irregular waves)
    // with A(w), B_rad(w) along the joint from pto_hydro.dat and the effective rigid-body mass
    // at the joint m = 1/(J^T M^-1 J); the damper (X 502) gets B_opt, or one is added
    double T = p->X506_T;

    if(T<=0.0)
    T = (p->B93==1) ? p->B93_2 : p->wT;

    pto_hydro_table tab;

    if(T<=0.0 || p->X506_C<=0.0 || !tab.read("pto_hydro.dat"))
    {
        if(p->mpirank==0)
        cout<<"6DOF PTO: X 506 needs T > 0 (or B 93), C > 0 and pto_hydro.dat (w A B_rad) in the case folder"<<endl;

        pgc->final();
        exit(1);
    }

    const double w = 2.0*PI/T;

    Eigen::Matrix<double,6,6> M = Eigen::Matrix<double,6,6>::Zero();
    M.block<3,3>(0,0) = Mass_fb*Eigen::Matrix3d::Identity();
    M.block<3,3>(3,3) = R_*I_*R_.transpose();

    const Eigen::Matrix<double,6,1> &J = pto_joint_.J;
    const double m = 1.0/J.dot(M.partialPivLu().solve(J));

    const double A = tab.A(w), B = tab.B(w);
    const double Bopt = pto_bopt(w, m, A, B, p->X506_C);

    if(!pto_.set_param("damper","B",Bopt))
    {
    pto_.elem.insert(pto_.elem.begin(), std::make_unique<pto_damper>(Bopt));     // first, as with X 502
    }

    if(p->mpirank==0)
    {
    cout<<"6DOF PTO tuned passive: T "<<T<<" s  m "<<m<<" kg  A "<<A<<" kg  B_rad "<<B<<" kg/s  C "<<p->X506_C<<" N/m  ->  B_opt "<<Bopt<<" N s/m"<<endl;

    if(w<tab.wmin() || w>tab.wmax())
    cout<<"6DOF PTO tuned passive: w = "<<w<<" outside the table, end values used"<<endl;
    }
}
