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

#ifndef SIXDOF_BRAIN_H_
#define SIXDOF_BRAIN_H_

#include<vector>
#include<Eigen/Dense>

//  The brain of a deforming 6DOF body (X 193): closes the loop between the senses of the fish and the
//  gait control (X 199). Once per half tail beat (the decision times) it reads what the body sensed
//  since the last decision and sets the targets of frequency, amplitude factor and curvature bias.
//
//  Senses (means over the time since the last decision, weighted with dt):
//   - the body axis and the lateral axis in the plane of the undulation (3D: x-y, 2D: x-z), the speed
//     along the axis, the position of the CoG (vision of the surroundings, over ground)
//   - the lateral line (X 194): pressure on both flanks per station
//   - vision (D > 0): 7 rays from the head in the plane of the undulation (0, +-25, +-50, +-75 deg to
//     the axis), the distance to the walls of the domain and to the solids (their level set) up to D,
//     at the decision times
//
//  The current and the stride: the fish knows its gait and sees its track over the ground. The speed of
//  the gait U_g follows f a L with the time 0.8 s; the ground velocity is s U_g ex + c. A Kalman filter
//  on the stride s (prior 0.5 L per beat) and the current c in the plane estimates both from the CoG at
//  the decisions (validation 11, 2D swimmer: s = 0.46, the current 0.17-0.26 m/s for the imposed 0.2 m/s,
//  in still water at most 0.08 m/s, mean 0.04 m/s). U_w = s U_g is the speed through the water. While the
//  fish glides (amplitude factor below 0.95 now or at the last decision) the estimate is not updated, and
//  U_g follows from the track.
//
//  Behaviours, by priority:
//   1. escape (th_esc > 0): a strong pressure signal at the head, |p+ - p-| or the excess of the head
//      pressure over the body mean, above th_esc q (q = rho U_w^2/2, through the water): a C-start away
//      from the side of the higher pressure (a 0.5, b -+2.5 in 0.08 s, after 0.15 s back to b = 0 at
//      1.5 times the frequency), then 2 s without a new escape. The signal of the own turn,
//      -0.5 b q at the head (validation 11), is taken off p+ - p- first. The flank lateral line senses
//      a wall only within ~0.1 L, so th_esc ~ 1.5: a startle at contact or a strong stimulus
//   2. avoidance (D > 0): the threat w = (1 - d/D)^2 of the closest ray; a turn to the side with the
//      longer rays (kept while the threat lasts), the bias blended from the target's to +-b_max with
//      min(1, w/0.25), the frequency times 1 - 0.3 min(1, w/0.25)
//   3. the target (X 192 id t x y z, linear in t): the wanted velocity over the ground is U_des towards
//      the target (slower within 1.5 L; beyond, U_des/cos of the angle between axis and target up to
//      60 deg, the closing speed on a detour) plus the velocity of the target as seen; the fish swims
//      w = wanted - c through the water. The heading: along w, the turn rate psi/tau (tau = 1.5 s,
//      psi the angle of w to the axis), the bias from omega = 1.3 b U_w/L. The frequency: |w|/(s a L),
//      no integrator; with w more than 1 rad off the axis at least the speed for the turn in 2.5 s at
//      b_max. Within 0.25 L: at the target. In a current the fish faces it and holds station.
//   4. glide: at the target (0.25 L), needing less than its slowest swimming speed through the water
//      (0.3 f0 s L), in a current slower than that and seeing nothing within D, the fish stops beating its tail and coasts (amplitude
//      factor 0, no bias); it beats again beyond 0.5 L of the target or when it needs 1.5 times that
//      speed: burst and coast at a station in still water or in a slow current
//  Without a target the fish keeps its frequency and swims straight.
//
//  The decisions: every half tail beat, at most 0.5 s apart. The brain is a function of the senses at
//  the decision times; its state and decisions go into the checkpoint, so a restart continues it exactly.

class sixdof_brain
{
public:

    struct target_point {double t, x, y, z;};

    sixdof_brain(double U_des, double b_max, double D_vis, double th_esc, double f0, double L, double rho,
                 double t_start, int lateral, const std::vector<target_point>&);

    // the senses of one time step: dt, the body axis and the lateral axis (inertial frame), the CoG,
    // the CoG velocity; the lateral line (per probe: s/L, side, pressure) when the body has one
    void sense(double dt, const Eigen::Vector3d &ax, const Eigen::Vector3d &lat, const Eigen::Vector3d &c,
               const Eigen::Vector3d &U, const std::vector<double> *ll_p, const std::vector<double> &ll_s,
               const std::vector<int> &ll_side);

    // a decision at time t (current gait: frequency f, amplitude factor a, bias b): new targets
    struct decision {bool on = false; double t, f, a, b, ts; bool second = false; double t2, f2, a2, b2, ts2;};
    bool due(double t) const {return t>=t_start_ && t>=next_;}
    decision decide(double t, double f, double a, double b, const std::vector<double> *vis);
    
    // vision: its range [m] (0: none) and the angles of the rays to the body axis [rad]
    double vision_range() const {return Dvis_*L_;}
    static const std::vector<double> &rays();

    // the log line of the last decision
    struct record {double t, tx, ty, tz, dist, psi, U, Ff, Fa, Fb, dpl, phd, dvis, s, cx, cy, cz, Uw; int mode;};
    const record &last() const {return rec_;}
    
    // what the brain perceived at the last decision (for an external brain, X 191): the mean axes in
    // the plane of the undulation, the speed, the CoG, the target, the mean pressures of the lateral
    // line, the vision rays, the estimate of the stride
    struct perception
    {
        Eigen::Vector3d ex, el, c, T;
        double U = 0.0, stride = 0.0;
        bool target = false;
        std::vector<double> llp, rays;
    };
    const perception &perceived() const {return perc_;}
    const std::vector<double> &ll_s() const {return lls_;}
    const std::vector<int> &ll_side() const {return llside_;}
    
    // an external brain (X 191) replaces the decision: the log record and the next decision time
    void external(const decision&, double next);

    // checkpoint
    void ckpt_state(std::vector<double>&) const;
    bool ckpt_restore(const std::vector<double>&);

    double fmax() const {return 2.0*f0_;}

private:

    Eigen::Vector3d target(double t) const;

    const double Udes_, bmax_, Dvis_, the_, f0_, L_, rho_, tstart_;
    const int lat_;
    std::vector<target_point> tp_;

    // accumulated senses since the last decision
    double T_ = 0.0, Uax_ = 0.0;
    Eigen::Vector3d ax_ = Eigen::Vector3d::Zero(), lt_ = Eigen::Vector3d::Zero(), c_ = Eigen::Vector3d::Zero();
    std::vector<double> llp_;
    double llT_ = 0.0;

    // state
    double t_start_ = 0.0, next_ = 0.0, refr_ = -1.0, side_ = 0.0;

    // the estimate of the stride and the current (Kalman filter on s, c in the plane), the speed of the
    // gait, the time, CoG and bias of the last decision, the target seen and its velocity
    void kf_update(const double H[2][3], const double y[2], double R, double dt);
    bool kf_on_ = false, Tset_ = false, glide_ = false;
    double kf_x_[3] = {0.0, 0.0, 0.0}, kf_P_[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    double ug_ = 0.0, tdec_ = 0.0, bdec_ = 0.0, adec_ = 1.0;
    Eigen::Vector3d pdec_ = Eigen::Vector3d::Zero(), Tprev_ = Eigen::Vector3d::Zero(), vT_ = Eigen::Vector3d::Zero();
    record rec_ = {};
    perception perc_;

    std::vector<double> lls_;
    std::vector<int> llside_;
};

#endif
