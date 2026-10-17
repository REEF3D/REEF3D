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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SIXDOF_EXTERNAL_H_
#define SIXDOF_EXTERNAL_H_

#include<vector>
#include<Eigen/Dense>

//  External rigid-body dynamics of the 6DOF bodies (REEF3D::CORAL, coral/reef3d, coral.txt fluid_body).
//
//  The fluid solver keeps the geometry, the body boundary conditions and the integration of the loads
//  on the hull; the external solver owns the equations of motion (joints, PTO, mooring) and returns the
//  pose and velocity of the bodies at the times the fluid solver needs them. Per fluid step t0 -> t1
//  (FNPF, fnpf_6DOF_external.cpp):
//
//   first stage:  fluid_loads(t0, F, A, false)      loads of the state at t0 (shared by both solvers)
//                 without predictor: advance(t0, t1), then state() for the stage poses
//                 with predictor:    predict() for the pose of the next stage (no step yet)
//   next stage:   with predictor:    fluid_loads(t_s, F, {}, true) (provisional, at the predicted pose),
//                                    then advance(t0, t1)
//   later stages: state(n, t) inside [t0, t1]
//
//  NHFLOW (6DOF_nhflow_external.cpp) integrates the loads at the end of each stage instead: the loads of
//  the initialisation are the first committed sample; per step the first stage takes predict() (or,
//  without predictor, advance(t0, t1) and state()), the loads at its end go over as provisional at
//  t0 + dt, then advance(t0, t1); the later stages take state(); the loads at the end of the step are
//  committed at t1. NHFLOW has no instantaneous added mass (A empty), its loads contain the reaction of
//  the body accelerations (CORAL: artificial added mass, fluid_body ... added_mass).
//
//  Loads: per body force and moment about the CoG, inertial frame, without gravity (the external solver
//  has its own); A: instantaneous added mass, block (n, m) = load on body
//  n from the accelerations of body m, index n*N + m (empty: none, or unchanged).
//  Bodies are numbered as the solver's 6DOF objects (0 .. N-1). Every rank makes the same calls with the
//  same (reduced) loads and gets the same states; the external solver needs no communication.
//  FNPF hands over the loads without the acceleration part and with A; NHFLOW with the acceleration part
//  and without A.
//
//  One-way bodies (one_way): the fluid solver keeps moving them itself (fixed, prescribed motion or its
//  own dynamics) and only hands over their loads, once per fluid step (FNPF: first stage, at t0; NHFLOW:
//  after the last stage, at t1, and the loads of the initialisation); the external solver does not
//  return a state for them.
//
//  The external solver registers itself (attach); the bodies look it up at their first step.

class sixdof_external
{
public:
    typedef Eigen::Matrix<double,6,1> vec6;
    typedef Eigen::Matrix<double,6,6> mat6;
    
    virtual ~sixdof_external() {}
    
    static sixdof_external* instance();
    static void attach(sixdof_external*);
    
    // body n is driven by the external solver (two-way), or only its loads go there (one_way)
    virtual bool drives(int n) const = 0;
    virtual bool one_way(int) const {return false;}
    // once, at the first step: mass, CoG, rotation (body -> inertial) and inertia (about the CoG, body
    // frame) of body n in REEF3D
    virtual void connect(int n, double mass, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Matrix3d &I) = 0;
    // the stage after the first enters the step (provisional loads)
    virtual bool predictor() const = 0;
    
    virtual void fluid_loads(double t, const std::vector<vec6> &F, const std::vector<mat6> &A, bool provisional) = 0;
    
    // distributed loads of a one-way body n: its hull forces binned along the axis a + s e, 0 <= s <= L,
    // in ns strips (false: none); the strip forces (3 per strip, this rank's part of the hull, called on
    // all ranks) at time t, before fluid_loads of the same time
    virtual bool strips(int, Eigen::Vector3d&, Eigen::Vector3d&, double&, int&) const {return false;}
    virtual void fluid_strips(int, double, const std::vector<double>&) {}
    virtual void advance(double t0, double t1) = 0;
    
    // pose and velocity of body n at time t (CoG, rotation body -> inertial, velocity, angular velocity
    // in the inertial frame): predicted from the current state (predict) or from the last step (state)
    virtual void predict(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) = 0;
    virtual void state(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) = 0;
};

#endif
