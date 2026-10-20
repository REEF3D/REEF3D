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
// REEF3D::CORAL: mooring lines between two attachment points (anchor on the ground: body -1, fairlead
// on a body; or between two bodies).
//
// mooring_qs: quasi-static elastic catenary (coral_catenary) in the vertical plane through both ends; the
// force follows the current positions instantly. The line has a seabed portion (friction CB) only if its
// lower end is a ground end on the seabed; otherwise (lines between two bodies, an anchor above the
// seabed, CB < 0) it hangs free between its ends.
//
// The seabed of a line (default_seabed): the height of its ground end (the lower one if both are on the
// ground); a line between two bodies (a shared line, a bridle) has none.
//
// mooring_line: dynamic lumped-mass line (MoorDyn type): N segments, the N-1 inner nodes carry mass,
// added mass (normal and tangential), submerged weight, axial stiffness EA (no compression) and axial
// damping BA (N s), drag normal and tangential against the water velocity (waves and current),
// Froude-Krylov inertia of the water, seabed contact (penalty kb [Pa/m], damping cb [Pa s/m],
// Coulomb friction mu with the regularization velocity v_eps). The inner nodes are advanced with an
// explicit Runge-Kutta scheme (RK4) and sub-steps inside the step of the multibody system, the end
// points moving linearly from the converged to the trial positions; eval() integrates from the
// converged line state (cached per trial state), commit() keeps the converged line. The force on an
// attached body is the tension of the end segment plus the weight and drag of the end half-segment.
// The line starts from the catenary profile (no seabed: straight).
//
// Wave kinematics (free surface, velocity, acceleration) are taken at the converged nodes at the start
// and at the end of the multibody step and are linear in time in between (the nodes move little against
// a wave length within a step); exact_kinematics evaluates them at every node position of every
// Runge-Kutta stage instead (costly with many wave components). settle() relaxes the line in the water
// at rest (the surface of the start, no wave velocities; the current acts).

#ifndef CORAL_MOORING_H_
#define CORAL_MOORING_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_state.h"
#include"coral/mooring/coral_catenary.h"
#include"coral/hydro/coral_waves.h"
#include<limits>
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    struct mooring_end
    {
        int body = -1;               // -1: ground
        Vec3 p = Vec3::Zero();       // global position at creation
    };

    // the seabed height of a line between a and b: the z of its ground end (the lower of two), none
    // (-1e30) between two bodies
    double default_seabed(const mooring_end &a, const mooring_end &b);

    struct mooring_props
    {
        double L = 0.0;              // unstretched length [m]
        double d = 0.1;              // volume-equivalent diameter [m]
        double mass = 0.0;           // mass per length in air [kg/m]
        double EA = 0.0;             // axial stiffness [N]
        double BA = 0.0;             // axial damping [N s]
        double CB = 0.0;             // seabed friction of the quasi-static line (< 0: no seabed)
        double Cdn = 1.2, Cdt = 0.05, Can = 1.0, Cat = 0.0;
        double kb = 3.0e6, cb = 3.0e5;     // seabed stiffness [Pa/m], damping [Pa s/m]
        double mu = 0.0, v_eps = 0.01;     // seabed friction (dynamic line)
        double rho = 1025.0;
    };

    class mooring_qs : public load
    {
    public:
        // z_seabed: NaN for default_seabed
        mooring_qs(const model&, const mooring_end &a, const mooring_end &b, const mooring_props&, const std::string &name="mooring",
                   double z_seabed=std::numeric_limits<double>::quiet_NaN());
        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(ba_); b.push_back(bb_);}
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // force on end b (global) and on end a for the current positions
        void forces(const model&, Vec3 &Fa, Vec3 &Fb) const;
        // restart from a checkpoint (coral_state.h)
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;

    private:
        std::string name_;
        int ba_, bb_;
        Vec3 sa_, sb_;
        mooring_props pr_;
        catenary_line cl_;
        int bed_end_ = -1;             // 0: a, 1: b is the ground end on the seabed (seabed portion), -1: none
        mutable double H_ = -1.0, V_ = 0.0;
    };

    class mooring_line : public load
    {
    public:
        mooring_line(const model&, const mooring_end &a, const mooring_end &b, const mooring_props&, int nseg,
                     std::shared_ptr<const wave_kinematics> waves=nullptr, const Vec3 &current=Vec3::Zero(),
                     double z_seabed=-1.0e30, double dt_line=0.0, const std::string &name="line");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(ba_); b.push_back(bb_);}
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        void kinematics_points(const model&, const void *source, std::vector<Vec3>&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        int nseg() const {return nseg_;}
        double dt_line() const {return dt_;}
        const std::vector<Vec3> &nodes() const {return r_;}     // converged node positions (0 .. N)
        double tension(int seg) const;                          // converged segment tension
        // end forces of the converged line (on the attached bodies)
        void end_forces(Vec3 &Fa, Vec3 &Fb) const {Fa = Fa_; Fb = Fb_;}

        // relax the line with the ends fixed (static shape), extra damping, before the run, in the water at
        // rest
        void settle(const model&, double t_settle);
        // wave kinematics at every Runge-Kutta stage (true) or once per step at the converged nodes (default)
        void set_exact_kinematics(bool b) {exact_kin_ = b;}
        // restart from a checkpoint (coral_state.h)
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;

    private:
        struct state {std::vector<Vec3> r, v;};
        // wave kinematics at the converged nodes at time t (gen: the converged line they belong to)
        struct node_kin {double t = 0.0; long gen = -1; std::vector<double> eta; std::vector<Vec3> u, du;};
        const node_kin &kinematics(int slot, double t) const;
        // the water at node i, position x, time t, for the integration from t0 to t1
        void water(int i, const Vec3 &x, double t, double &eta, Vec3 &u, Vec3 &du) const;
        void accel(const state&, const Vec3 &ra, const Vec3 &rb, const Vec3 &va, const Vec3 &vb, double t,
                   std::vector<Vec3> &a, Vec3 &Fa, Vec3 &Fb, double extra_damp) const;
        // from the converged line at t0 to t1; at_rest: settling in the water at rest
        void integrate(state&, const Vec3 &ra0, const Vec3 &ra1, const Vec3 &rb0, const Vec3 &rb1,
                       double t0, double t1, Vec3 &Fa, Vec3 &Fb, double extra_damp, bool at_rest=false) const;

        std::string name_;
        int ba_, bb_;
        Vec3 sa_, sb_;
        mooring_props pr_;
        int nseg_;
        double l_;                     // unstretched segment length
        double dt_;
        std::shared_ptr<const wave_kinematics> waves_;
        Vec3 current_;
        double zbed_;
        Vec3 g_;

        // converged state: inner nodes (and end positions) at t_
        std::vector<Vec3> r_, v_;
        Vec3 ra_, rb_;
        double t_ = 0.0;
        Vec3 Fa_ = Vec3::Zero(), Fb_ = Vec3::Zero();
        long gen_ = 0;                 // changes with the converged line

        // wave kinematics: per step at the converged nodes (slots: start, end of the step) or exact
        bool exact_kin_ = false;
        mutable node_kin kin_[2];
        mutable int kmode_ = 0;        // 0: no waves, 1: exact, 2: linear between kin_[0] and kin_[1], 3: at rest (kin_[0] surface)
        mutable double kt0_ = 0.0, kt1_ = 0.0;

        // work arrays
        mutable std::vector<Vec3> T_;
        mutable std::vector<Vec3> w_[4];
        mutable state ws_;

        // trial cache
        mutable Vec3 cra_, crb_;
        mutable double ct_ = -1.0e300;
        mutable Vec3 cFa_, cFb_;
        mutable bool cvalid_ = false;
    };
}

#endif
