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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_rigidbody.h"
#include"definitions.h"
#include<cmath>

sixdof_rigidbody::sixdof_rigidbody() : twoD(false), mass(1.0), phi(0.0), theta(0.0), psi(0.0),
                                       dtn1(0.0), dtn2(0.0), dtn3(0.0)
{
    for(int n=0; n<6; ++n)
    dof[n] = 1;

    for(int n=0; n<3; ++n)
    Cdamp_t[n] = Cdamp_r[n] = 0.0;

    I.setIdentity();

    reset();

    pk.setZero(); dpk.setZero();
    ck.setZero(); dck.setZero();
    hk.setZero(); dhk.setZero();
    ek.setZero(); dek.setZero();
    pn1.setZero(); pn2.setZero(); pn3.setZero();
    cn1.setZero(); cn2.setZero(); cn3.setZero();
    hn1.setZero(); hn2.setZero(); hn3.setZero();
    en1.setZero(); en2.setZero(); en3.setZero();
    E.setZero(); G.setZero(); Gdot.setZero();
    Rinv.setZero();

    for(int s=0; s<3; ++s)
    {
    rk4_p[s].setZero();
    rk4_c[s].setZero();
    rk4_h[s].setZero();
    rk4_e[s].setZero();
    }
}

// ---------------------------------------------------------------------------------------------
// set-up
// ---------------------------------------------------------------------------------------------

void sixdof_rigidbody::reset()
{
    R << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    e << 0.0, 0.0, 0.0, 0.0;
    p << 0.0, 0.0, 0.0;
    c << 0.0, 0.0, 0.0;
    h << 0.0, 0.0, 0.0;

    dp   << 0.0, 0.0, 0.0;
    dpn1 << 0.0, 0.0, 0.0;
    dpn2 << 0.0, 0.0, 0.0;
    dpn3 << 0.0, 0.0, 0.0;

    dc   << 0.0, 0.0, 0.0;
    dcn1 << 0.0, 0.0, 0.0;
    dcn2 << 0.0, 0.0, 0.0;
    dcn3 << 0.0, 0.0, 0.0;

    dh   << 0.0, 0.0, 0.0;
    dhn1 << 0.0, 0.0, 0.0;
    dhn2 << 0.0, 0.0, 0.0;
    dhn3 << 0.0, 0.0, 0.0;

    de   << 0.0, 0.0, 0.0, 0.0;
    den1 << 0.0, 0.0, 0.0, 0.0;
    den2 << 0.0, 0.0, 0.0, 0.0;
    den3 << 0.0, 0.0, 0.0, 0.0;

    omega_B << 0.0, 0.0, 0.0;
    omega_I << 0.0, 0.0, 0.0;

    phi = theta = psi = 0.0;

    F << 0.0, 0.0, 0.0;
    M << 0.0, 0.0, 0.0;
}

void sixdof_rigidbody::quaternion_from_euler()
{
	// Goldstein p. 604
	e(0) =
		 cos(0.5*phi)*cos(0.5*theta)*cos(0.5*psi)
		+ sin(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e(1) =
		 sin(0.5*phi)*cos(0.5*theta)*cos(0.5*psi)
		- cos(0.5*phi)*sin(0.5*theta)*sin(0.5*psi);
	e(2) =
		 cos(0.5*phi)*sin(0.5*theta)*cos(0.5*psi)
		+ sin(0.5*phi)*cos(0.5*theta)*sin(0.5*psi);
	e(3) =
		 cos(0.5*phi)*cos(0.5*theta)*sin(0.5*psi)
		- sin(0.5*phi)*sin(0.5*theta)*cos(0.5*psi);
}

void sixdof_rigidbody::init_history()
{
    en1 = e;
    en2 = e;
    en3 = e;
    ek  = e;

    cn1 = c;
    cn2 = c;
    cn3 = c;
    ck  = c;

    pn1 = p;
    pn2 = p;
    pn3 = p;
    pk  = p;

    hn1 = h;
    hn2 = h;
    hn3 = h;
    hk  = h;

    dpk = dp;
    dck = dc;
    dhk = dh;
    dek = de;
}

bool sixdof_rigidbody::fixed(int n) const
{
    if(twoD && (n==1 || n==3 || n==5))
    return true;

    return (dof[n]!=1);
}

// ---------------------------------------------------------------------------------------------
// kinematics
// ---------------------------------------------------------------------------------------------

void sixdof_rigidbody::quat_matrices()
{
    // Shivarama PhD thesis, p. 19
    E << -e(1), e(0), -e(3), e(2),
         -e(2), e(3), e(0), -e(1),
         -e(3), -e(2), e(1), e(0);

    G << -e(1), e(0), e(3), -e(2),
         -e(2), -e(3), e(0), e(1),
         -e(3), e(2), -e(1), e(0);

    R = E*G.transpose();
    Rinv = R.inverse();
}

void sixdof_rigidbody::euler_angles()
{
	// around z-axis
	psi = atan2(2.0*(e(1)*e(2) + e(3)*e(0)), 1.0 - 2.0*(e(2)*e(2) + e(3)*e(3)));

	// around new y-axis
	double arg = 2.0*(e(0)*e(2) - e(1)*e(3));

	if(fabs(arg) >= 1.0)
	theta = (arg>=0.0 ? 1.0 : -1.0)*PI/2.0;

	else
	theta = asin(arg);

	// around new x-axis
	phi = atan2(2.0*(e(2)*e(3) + e(1)*e(0)), 1.0 - 2.0*(e(1)*e(1) + e(2)*e(2)));
}

void sixdof_rigidbody::update_omega()
{
    omega_B = I.inverse()*h;
    omega_I = R*omega_B;
}

void sixdof_rigidbody::velocity(Eigen::Matrix<double, 6, 1> &u) const
{
    // translation: momentum for free DOFs, the prescribed velocity for prescribed DOFs
    if(dof[0]==0)
    u(0) = 0.0;
    if(dof[0]==1)
    u(0) = p(0)/mass;
    if(dof[0]==2)
    u(0) = dc(0);

    if(dof[1]==0 || twoD)
    u(1) = 0.0;
    if(dof[1]==1 && !twoD)
    u(1) = p(1)/mass;
    if(dof[1]==2)
    u(1) = dc(1);

    if(dof[2]==0)
    u(2) = 0.0;
    if(dof[2]==1)
    u(2) = p(2)/mass;
    if(dof[2]==2)
    u(2) = dc(2);

    // rotation
    if(twoD)
    {
    u(3) = 0.0;
    u(4) = omega_I(1);
    u(5) = 0.0;
    }

    if(!twoD)
    {
    u(3) = omega_I(0);
    u(4) = omega_I(1);
    u(5) = omega_I(2);
    }
}

// ---------------------------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------------------------

void sixdof_rigidbody::assemble_loads(const double *Fext)
{
    // Fext: X, Y, Z, K, M, N (inertial frame, moments about the CoG) without damping
    F << 0.0, 0.0, 0.0;
    M << 0.0, 0.0, 0.0;

    for(int n=0; n<3; ++n)
    if(dof[n]==1)
    F(n) = Fext[n] - Cdamp_t[n]*p(n)/mass;

    for(int n=0; n<3; ++n)
    if(dof[n+3]==1)
    M(n) = Fext[n+3] - Cdamp_r[n]*omega_I(n);
}

// ---------------------------------------------------------------------------------------------
// right-hand side
// ---------------------------------------------------------------------------------------------

void sixdof_rigidbody::derivatives_trans()
{
    dp = F;       // d(linear momentum)/dt = force
    dc = p/mass;  // d(CoG position)/dt = velocity
}

void sixdof_rigidbody::derivatives_rot()
{
    quat_matrices();

    // RHS of e
    de = 0.5*G.transpose()*I.inverse()*h;

    // RHS of h: moment transformed into the body-fixed system (Shivarama and Schwab)
    Gdot << -de(1), de(0), de(3),-de(2),
            -de(2),-de(3), de(0), de(1),
            -de(3), de(2),-de(1), de(0);

    dh = 2.0*Gdot*G.transpose()*h + Rinv*M;
}

// ---------------------------------------------------------------------------------------------
// stage updates
// ---------------------------------------------------------------------------------------------

void sixdof_rigidbody::stage_rk2(int iter, double dt)
{
    if(iter==0)
    {
        pk = p;
        ck = c;
        hk = h;
        ek = e;

        p = pk + dt*dp;
        c = ck + dt*dc;
        h = hk + dt*dh;
        e = ek + dt*de;
        e.normalize();
    }

    if(iter==1)
    {
        p = 0.5*pk + 0.5*p + 0.5*dt*dp;
        c = 0.5*ck + 0.5*c + 0.5*dt*dc;
        h = 0.5*hk + 0.5*h + 0.5*dt*dh;
        e = 0.5*ek + 0.5*e + 0.5*dt*de;
        e.normalize();
    }
}

void sixdof_rigidbody::stage_rk3(int iter, double dt)
{
    if(iter==0)
    {
        pk = p;
        ck = c;
        hk = h;
        ek = e;

        p = pk + dt*dp;
        c = ck + dt*dc;
        h = hk + dt*dh;
        e = ek + dt*de;
        e.normalize();
    }

    if(iter==1)
    {
        p = 0.75*pk + 0.25*p + 0.25*dt*dp;
        c = 0.75*ck + 0.25*c + 0.25*dt*dc;
        h = 0.75*hk + 0.25*h + 0.25*dt*dh;
        e = 0.75*ek + 0.25*e + 0.25*dt*de;
        e.normalize();
    }

    if(iter==2)
    {
        p = (1.0/3.0)*pk + (2.0/3.0)*p + (2.0/3.0)*dt*dp;
        c = (1.0/3.0)*ck + (2.0/3.0)*c + (2.0/3.0)*dt*dc;
        h = (1.0/3.0)*hk + (2.0/3.0)*h + (2.0/3.0)*dt*dh;
        e = (1.0/3.0)*ek + (2.0/3.0)*e + (2.0/3.0)*dt*de;
        e.normalize();
    }
}

void sixdof_rigidbody::stage_rkls3(double gamma, double zeta, double dt)
{
    p = p + gamma*dt*dp + zeta*dt*dpk;
    c = c + gamma*dt*dc + zeta*dt*dck;
    h = h + gamma*dt*dh + zeta*dt*dhk;
    e = e + gamma*dt*de + zeta*dt*dek;
    e.normalize();

    dpk = dp;
    dck = dc;
    dhk = dh;
    dek = de;
}

void sixdof_rigidbody::stage_rk4(int iter, double dt)
{
    // classical RK4:
    // y_s = y_n + c_s*dt*K_s (c_s = 1/2, 1/2, 1),  y_n+1 = y_n + dt/6*(K1 + 2K2 + 2K3 + K4)
    if(iter==0)
    {
        pk = p;
        ck = c;
        hk = h;
        ek = e;
    }

    if(iter<3)
    {
        rk4_p[iter] = dp;
        rk4_c[iter] = dc;
        rk4_h[iter] = dh;
        rk4_e[iter] = de;

        const double cs = (iter==2) ? 1.0 : 0.5;

        p = pk + cs*dt*dp;
        c = ck + cs*dt*dc;
        h = hk + cs*dt*dh;
        e = ek + cs*dt*de;
        e.normalize();
    }

    if(iter==3)
    {
        const double w = dt/6.0;

        p = pk + w*(rk4_p[0] + 2.0*rk4_p[1] + 2.0*rk4_p[2] + dp);
        c = ck + w*(rk4_c[0] + 2.0*rk4_c[1] + 2.0*rk4_c[2] + dc);
        h = hk + w*(rk4_h[0] + 2.0*rk4_h[1] + 2.0*rk4_h[2] + dh);
        e = ek + w*(rk4_e[0] + 2.0*rk4_e[1] + 2.0*rk4_e[2] + de);
        e.normalize();
    }
}

void sixdof_rigidbody::step_onestep(double dt)
{
    // Heun form with the derivatives of the step start: for prescribed motion (constant
    // derivatives) this is exact, for general loads it reduces to explicit Euler
    pk = p;
    ck = c;
    hk = h;
    ek = e;

    p = p + dt*dp;
    c = c + dt*dc;
    h = h + dt*dh;
    e = e + dt*de;
    e.normalize();

    p = 0.5*pk + 0.5*p + 0.5*dt*dp;
    c = 0.5*ck + 0.5*c + 0.5*dt*dc;
    h = 0.5*hk + 0.5*h + 0.5*dt*dh;
    e = 0.5*ek + 0.5*e + 0.5*dt*de;
    e.normalize();
}

void sixdof_rigidbody::save_history(double dtstage)
{
    dpn3 = dpn2;
    dpn2 = dpn1;
    dpn1 = dp;

    dcn3 = dcn2;
    dcn2 = dcn1;
    dcn1 = dc;

    dhn3 = dhn2;
    dhn2 = dhn1;
    dhn1 = dh;

    den3 = den2;
    den2 = den1;
    den1 = de;

    pn3 = pn2;
    pn2 = pn1;
    pn1 = p;

    cn3 = cn2;
    cn2 = cn1;
    cn1 = c;

    hn3 = hn2;
    hn2 = hn1;
    hn1 = h;

    en3 = en2;
    en2 = en1;
    en1 = e;

    dtn3 = dtn2;
    dtn2 = dtn1;
    dtn1 = dtstage;
}
