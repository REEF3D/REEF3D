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

#ifndef SIXDOF_RIGIDBODY_H_
#define SIXDOF_RIGIDBODY_H_

#include<Eigen/Dense>

//  Rigid-body core of the 6DOF module: state, kinematics and time integration of one body.
//
//  Solver independent: no lexer, fdm or ghostcell. The hydrodynamic model, mooring, PTO etc.
//  deliver the resultant load (F, M); prescribed motions (sixdof_motionext) act on the
//  derivatives between derivatives_*() and the stage update.
//
//  State (Shivarama & Schwab, quaternion formulation):
//   c  position of the centre of gravity, inertial frame
//   p  linear momentum, inertial frame
//   h  angular momentum, body-fixed frame
//   e  unit quaternion (Euler parameters), body -> inertial
//
//  Loads: F force in the inertial frame, M moment about the CoG in the inertial frame.
//
//  DOF modes (dof[0..5] = surge, sway, heave, roll, pitch, yaw):
//   0 fixed, 1 free, 2 prescribed (the velocity comes from sixdof_motionext)
//
//  A time step with an s-stage scheme is:
//   for each stage: (loads) -> assemble_loads() -> derivatives_trans() / derivatives_rot()
//                   -> (prescribed motion) -> stage_*() -> quat_matrices()
//   at the end:     save_history()

class sixdof_rigidbody
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    sixdof_rigidbody();

    // ---- configuration
    int dof[6];          // DOF modes, see above
    bool twoD;           // 2D (x-z plane): sway, roll and yaw do not exist
    double Cdamp_t[3];   // linear damping of the translations   [N s/m]   (X 26)
    double Cdamp_r[3];   // linear damping of the rotations      [N m s]   (X 25)

    // ---- mass properties
    double mass;
    Eigen::Matrix3d I;   // inertia tensor about the CoG, body-fixed frame

    // ---- state, stage start (k), history (n1..n3) and derivatives
    Eigen::Vector3d p, pk, pn1, pn2, pn3, dp, dpk, dpn1, dpn2, dpn3;
    Eigen::Vector3d c, ck, cn1, cn2, cn3, dc, dck, dcn1, dcn2, dcn3;
    Eigen::Vector3d h, hk, hn1, hn2, hn3, dh, dhk, dhn1, dhn2, dhn3;
    Eigen::Vector4d e, ek, en1, en2, en3, de, dek, den1, den2, den3;

    // ---- kinematics derived from e and h
    Eigen::Matrix<double, 3, 4> E, G, Gdot;
    Eigen::Matrix3d R, Rinv;   // rotation matrix body -> inertial and its inverse
    Eigen::Vector3d omega_B;   // angular velocity, body-fixed frame
    Eigen::Vector3d omega_I;   // angular velocity, inertial frame
    double phi, theta, psi;    // Euler angles (roll, pitch, yaw) [rad]

    // ---- resultant load on the body
    Eigen::Vector3d F, M;

    // ---- stage time steps of the last steps
    double dtn1, dtn2, dtn3;

    // ---- set-up
    void reset();                                       // zero state, derivatives, history and loads
    void quaternion_from_euler();                       // e from phi, theta, psi (Goldstein p. 604)
    void init_history();                                // stage and history copies = current state
    bool fixed(int) const;                              // DOF n not free (fixed or prescribed, or absent in 2D)

    // ---- kinematics
    void quat_matrices();                               // E, G, R, Rinv from e
    void euler_angles();                                // phi, theta, psi from e
    void update_omega();                                // omega_B, omega_I from h
    void velocity(Eigen::Matrix<double, 6, 1>&) const;  // (u, v, w, p, q, r) of the CoG, inertial frame

    // ---- loads
    // F, M from the external loads (inertial frame, about the CoG) and the linear damping;
    // fixed and prescribed DOFs get zero load
    void assemble_loads(const double*);

    // ---- right-hand side
    void derivatives_trans();                           // dp = F, dc = p/m
    void derivatives_rot();                             // de, dh (calls quat_matrices)

    // ---- stage updates
    void stage_rk2(int, double);                        // TVD RK2 (Heun), stage iter, dt
    void stage_rk3(int, double);                        // TVD RK3 (Shu-Osher)
    void stage_rkls3(double, double, double);           // low-storage RK3: gamma, zeta, dt
    void stage_rk4(int, double);                        // classical RK4, stages 0..3
    void step_onestep(double);                          // single step with frozen derivatives (prescribed motion)

    void save_history(double);                          // shift n1..n3, argument: stage time step

private:
    Eigen::Vector3d rk4_p[3], rk4_c[3], rk4_h[3];
    Eigen::Vector4d rk4_e[3];
};

#endif
