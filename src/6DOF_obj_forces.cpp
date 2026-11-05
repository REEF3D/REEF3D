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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"

void sixdof_obj_cfd::hydrodynamic_forces_cfd(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter, bool finalize)
{
    if(p->X60==1)
    forces_stl(p,a,pgc,uvel,vvel,wvel,iter,finalize);
    
    if(p->X60==2)
    forces_lsm(p,a,pgc,uvel,vvel,wvel,iter,finalize);
}

void sixdof_obj::update_forces(lexer *p)
{
    // Forces in inertial system: external loads, linear damping, DOF modes
    double Fext[6] = {Xext + Xe, Yext + Ye, Zext + Ze, Kext + Ke, Mext + Me, Next + Ne};
    
    // load models (ship module): evaluated with the state of the stage
    bool am_load = false;
    Eigen::Matrix<double,6,6> A_load;
    
    if(!pload.empty())
    {
        double Fl[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        
        for(size_t ql=0; ql<pload.size(); ++ql)
        pload[ql]->add_load(p,rb,geom,pfluid,Fl);
        
        // hydrodynamic loads of the coupling that a model replaces (e.g. MMG manoeuvring)
        double w[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
        
        for(size_t ql=0; ql<pload.size(); ++ql)
        pload[ql]->fluid_mask(w);
        
        const double Fh[6] = {Xe, Ye, Ze, Ke, Me, Ne};
        
        for(int qn=0; qn<6; ++qn)
        {
            if(w[qn]!=1.0)
            Fext[qn] -= (1.0 - w[qn])*Fh[qn];
            
            Fext[qn] += Fl[qn];
        }
        
        // added mass of the models (body frame) -> inertial frame
        Eigen::Matrix<double,6,6> Ab = Eigen::Matrix<double,6,6>::Zero();
        
        for(size_t ql=0; ql<pload.size(); ++ql)
        am_load = pload[ql]->added_mass(rb,Ab) || am_load;
        
        if(am_load)
        {
            Eigen::Matrix<double,6,6> T = Eigen::Matrix<double,6,6>::Zero();
            T.block<3,3>(0,0) = rb.R;
            T.block<3,3>(3,3) = rb.R;
            A_load = T*Ab*T.transpose();
        }
    }
    
    rb.assemble_loads(Fext);
    
    if(Ffb_(0)!=Ffb_(0))
    cout<<"Ffb_(0)....###"<<endl;
    
    if(Ffb_(1)!=Ffb_(1))
    cout<<"Ffb_(1)....###"<<endl;
    
    if(Ffb_(2)!=Ffb_(2))
    cout<<"Ffb_(2)....###"<<endl;
    
    
    if(Mfb_(0)!=Mfb_(0))
    cout<<"Mfb_(0)....###"<<endl;
    
    if(Mfb_(1)!=Mfb_(1))
    cout<<"Mfb_(1)....###"<<endl;
    
    if(Mfb_(2)!=Mfb_(2))
    cout<<"Mfb_(2)....###"<<endl;
    
    // FNPF: instantaneous added mass (and implicit PTO terms, X 500 2) on the left-hand side
    if(am_on_ || am_load || (pto_on_ && pto_implicit_))
    apply_added_mass(p, am_load ? &A_load : nullptr);
    else if(am_cfd_>0.0)
    cfd_added_mass(p);
}

// CFD: the loads of the fluid contain the reaction of the fluid to the acceleration of the body of
// the last stage (the pressure of its projection).  For a light body the explicit coupling then
// oscillates from stage to stage and grows (the loads of the finest grid of the mesh refinement
// with G 7 0, a box of density 500 in heave: validation 51).  An artificial added mass A on both
// sides, (M + A) a = F + A a_prev, damps it and is consistent where the acceleration is smooth
// (a = a_prev): A = am_cfd_ times the mass of the displaced fluid in translation, times the
// inertia of the body with the density of the fluid in rotation
void sixdof_obj::cfd_added_mass(lexer *p)
{
    const Eigen::Matrix3d II = R_*I_*R_.transpose();
    const Eigen::Vector3d w = omega_I;
    const Eigen::Vector3d gyro = w.cross(II*w);
    const double rf = p->W1;
    const double rb = (Vfb>0.0) ? Mass_fb/Vfb : rf;

    Eigen::Matrix<double,6,6> A = Eigen::Matrix<double,6,6>::Zero();
    
    if(am_ell_on_)
    {
    // the added mass of the ellipsoid and the fluid inside the hull, which the direct forcing
    // accelerates with the body (its reaction is in the loads as well)
    Eigen::Matrix<double,6,6> T = Eigen::Matrix<double,6,6>::Zero();
    T.block<3,3>(0,0) = R_;
    T.block<3,3>(3,3) = R_;
    A = T*am_ell_*T.transpose();
    A.block<3,3>(0,0) += rf*Vfb*Eigen::Matrix3d::Identity();
    A.block<3,3>(3,3) += (rf/rb)*II;
    A *= am_cfd_;
    }
    
    else
    {
    A.block<3,3>(0,0) = am_cfd_*rf*Vfb*Eigen::Matrix3d::Identity();
    A.block<3,3>(3,3) = am_cfd_*(rf/rb)*II;
    }

    Eigen::Matrix<double,6,6> L = A;
    L.block<3,3>(0,0) += Mass_fb*Eigen::Matrix3d::Identity();
    L.block<3,3>(3,3) += II;

    Eigen::Matrix<double,6,1> r;
    r.head<3>() = Ffb_;
    r.tail<3>() = Mfb_ - gyro;
    r += A*am_prev_;

    bool fixed[6];
    for(int n=0; n<6; ++n)
    {
        fixed[n] = p_fixed_dof(p,n);
        if(fixed[n])
        {
            L.row(n).setZero();
            L.col(n).setZero();
            L(n,n) = 1.0;
            r(n) = 0.0;
        }
    }

    const Eigen::Matrix<double,6,1> acc = L.partialPivLu().solve(r);
    am_prev_ = acc;
    finish_added_mass(acc,II,gyro,fixed);
}

void sixdof_obj::added_mass_ellipsoid(lexer *p)
{
    // half extents of the hull in the body frame (relative to the CoG): the semi-axes a, b, c of
    // the equivalent ellipsoid along x, y, z
    double lo[3] = {1.0e20,1.0e20,1.0e20}, hi[3] = {-1.0e20,-1.0e20,-1.0e20};
    
    for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        const double v[3] = {tri_x0[n][q], tri_y0[n][q], tri_z0[n][q]};
        
        for(int d=0; d<3; ++d)
        {
        lo[d] = MIN(lo[d],v[d]);
        hi[d] = MAX(hi[d],v[d]);
        }
    }
    
    const double a = 0.5*(hi[0]-lo[0]), b = 0.5*(hi[1]-lo[1]), c = 0.5*(hi[2]-lo[2]);
    const double rf = p->W1;
    
    // added mass of the ellipsoid about its centre (diagonal in the body frame)
    Eigen::Matrix<double,6,1> m = Eigen::Matrix<double,6,1>::Zero();
    
    if(p->j_dir==0)
    {
        // 2D (x-z plane): elliptic cylinder of the span B = 2b (Lamb): surge rho pi c^2 B, heave
        // rho pi a^2 B, pitch rho pi (a^2 - c^2)^2 B/8
        const double B = 2.0*b;
        m(0) = rf*PI*c*c*B;
        m(2) = rf*PI*a*a*B;
        m(4) = rf*PI*pow(a*a - c*c,2.0)*B/8.0;
    }
    
    else
    {
        // ellipsoid (Lamb 1932, 373): alpha0 = abc int_0^inf dl/((a^2+l) D), D = sqrt((a^2+l)(b^2+l)(c^2+l)),
        // beta0, gamma0 likewise; l = m^2 e^s, trapezoidal rule in s
        const double mx = MAX(a,MAX(b,c));
        double al = 0.0, be = 0.0, ga = 0.0;
        const int N = 8000;
        const double s0 = -40.0, s1 = 40.0, ds = (s1 - s0)/double(N);
        
        for(int q=0; q<=N; ++q)
        {
            const double l = mx*mx*exp(s0 + q*ds);
            const double D = sqrt((a*a + l)*(b*b + l)*(c*c + l));
            const double w = ((q==0 || q==N) ? 0.5 : 1.0)*ds*l/D;
            al += w/(a*a + l);
            be += w/(b*b + l);
            ga += w/(c*c + l);
        }
        
        al *= a*b*c;
        be *= a*b*c;
        ga *= a*b*c;
        
        const double V = (4.0/3.0)*PI*a*b*c;
        
        m(0) = rf*V*al/(2.0 - al);
        m(1) = rf*V*be/(2.0 - be);
        m(2) = rf*V*ga/(2.0 - ga);
        
        // rotations: (rho V/5) (b^2 - c^2)^2 (gamma0 - beta0)/(2 (b^2 - c^2) + (b^2 + c^2)(beta0 - gamma0)), cyclic
        auto rot = [&](double b2, double c2, double bb, double gg)
        {
            const double den = 2.0*(b2 - c2) + (b2 + c2)*(bb - gg);
            return (fabs(b2 - c2)>1.0e-12*(b2 + c2) && fabs(den)>1.0e-300) ? (rf*V/5.0)*pow(b2 - c2,2.0)*(gg - bb)/den : 0.0;
        };
        
        m(3) = rot(b*b,c*c,be,ga);
        m(4) = rot(c*c,a*a,ga,al);
        m(5) = rot(a*a,b*b,al,be);
    }
    
    // about the CoG: the velocity of the centre r of the ellipsoid is v + w x r = v - [r]x w, so
    // A = J^T Mt J + diag(0, Ic) with J = [1, -[r]x]
    const Eigen::Vector3d r(0.5*(hi[0]+lo[0]), 0.5*(hi[1]+lo[1]), 0.5*(hi[2]+lo[2]));
    Eigen::Matrix3d S;
    S <<  0.0, -r(2),  r(1),
         r(2),   0.0, -r(0),
        -r(1),  r(0),   0.0;
    
    Eigen::Matrix<double,3,6> J;
    J.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
    J.block<3,3>(0,3) = -S;
    
    const Eigen::Matrix3d Mt = m.head<3>().asDiagonal();
    am_ell_ = J.transpose()*Mt*J;
    am_ell_.block<3,3>(3,3) += Eigen::Matrix3d(m.tail<3>().asDiagonal());
    
    // 2D: sway, roll and yaw do not exist
    if(p->j_dir==0)
    for(int q=0; q<6; ++q)
    if(q==1 || q==3 || q==5)
    {
    am_ell_.row(q).setZero();
    am_ell_.col(q).setZero();
    }
    
    am_ell_on_ = true;
    
    if(p->mpirank==0)
    cout<<"6DOF body "<<n6DOF<<": stabilised coupling X 197 "<<p->X197<<", added mass of the equivalent "<<(p->j_dir==0 ? "elliptic cylinder" : "ellipsoid")
        <<" (semi-axes "<<a<<", "<<b<<", "<<c<<", centre "<<r.transpose()<<" from the CoG): "<<m(0)<<" "<<m(1)<<" "<<m(2)<<" kg, "<<m(3)<<" "<<m(4)<<" "<<m(5)<<" kg m^2 (mass "<<Mass_fb<<" kg)"<<endl;
}

void sixdof_obj::apply_added_mass(lexer *p, const Eigen::Matrix<double,6,6> *A_extra)
{
    // Rigid body with the instantaneous added mass A (inertial frame, moments about the CoG):
    //   [M I + A_tt   A_tr    ] [a    ]   [ F                ]
    //   [A_rt      I_I + A_rr ] [alpha] = [ Mo - w x (I_I w) ]
    // update_forces has already assembled F and Mo (hydrodynamic without the
    // acceleration part, gravity, mooring, damping). The kernel integrates
    // dp/dt = F and dh_B/dt = 2 Gdot G^T h + R^T Mo, so handing it
    //   F* = M a,   Mo* = I_I alpha + w x (I_I w)
    // gives dh_B/dt = I_B alpha_B and leaves get_trans/get_rot untouched.
    
    const Eigen::Matrix3d II = R_*I_*R_.transpose();
    const Eigen::Vector3d w = omega_I;
    const Eigen::Vector3d gyro = w.cross(II*w);
    
    // FNPF added mass (am_on_) and that of the load models (A_extra, e.g. MMG)
    Eigen::Matrix<double,6,6> L = Aadd_;
    
    if(A_extra!=nullptr)
    {
        if(!am_on_)
        L = *A_extra;
        else
        L += *A_extra;
    }
    
    L.block<3,3>(0,0) += Mass_fb*Eigen::Matrix3d::Identity();
    L.block<3,3>(3,3) += II;
    
    Eigen::Matrix<double,6,1> r;
    r.head<3>() = Ffb_;
    r.tail<3>() = Mfb_ - gyro;
    
    // PTO Jacobians (X 500 2)
    pto_implicit(p,L,r);

    bool fixed[6];
    
    for(int n=0; n<6; ++n)
    {
    fixed[n] = p_fixed_dof(p,n);
    
    if(fixed[n])
    {
    L.row(n).setZero();
    L.col(n).setZero();
    L(n,n) = 1.0;
    r(n) = 0.0;
    }
    }
    
    // several FNPF bodies: hand the block to the coupled solve (fnpf_6DOF::motion)
    if(am_coupled_)
    {
        Lc_ = L;
        rc_ = r;
        IIc_ = II;
        gyroc_ = gyro;

        for(int n=0; n<6; ++n)
        fixedc_[n] = fixed[n];

        am_assembled_ = true;
        return;
    }

    const Eigen::Matrix<double,6,1> acc = L.partialPivLu().solve(r);
    
    finish_added_mass(acc,II,gyro,fixed);
}

void sixdof_obj::finish_added_mass(const Eigen::Matrix<double,6,1> &acc, const Eigen::Matrix3d &II,
                                   const Eigen::Vector3d &gyro, const bool *fixed)
{
    Ffb_ = Mass_fb*acc.head<3>();
    Mfb_ = II*acc.tail<3>() + gyro;
    
    // fixed DOFs keep the kernel's convention of zero load
    for(int n=0; n<3; ++n)
    {
    if(fixed[n])
    Ffb_(n) = 0.0;
    
    if(fixed[n+3])
    Mfb_(n) = 0.0;
    }
}

bool sixdof_obj::p_fixed_dof(lexer *p, int n)
{
    // free DOF: X11 flag 1 (2 = prescribed via motionext, 0 = fixed)
    // 2D: sway, roll and yaw do not exist
    return rb.fixed(n);
}

