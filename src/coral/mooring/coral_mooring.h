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
// Froude-Krylov inertia of the water (taken by the submerged fraction of a node's piece of line, linear
// over its vertical extent, so nodes cross the free surface smoothly), seabed contact (penalty kb [Pa/m],
// damping cb [Pa s/m], Coulomb friction mu across and mu_axial along the line with the regularization
// velocity v_eps). The inner nodes are advanced with an
// explicit Runge-Kutta scheme (RK4) and sub-steps inside the step of the multibody system, the end
// points moving on cubic Hermite curves from the converged to the trial positions and velocities of the
// attachment points; eval() integrates from the converged line state (cached per trial state), commit()
// keeps the converged line. The force on an attached body is the tension of the end segment plus the
// weight, drag and Froude-Krylov force of the end half-segment, minus its inertia: mass and added mass
// (normal and tangential, submerged) times the acceleration of the attachment point (as the coupled
// points of MoorDyn), so the body carries the end node. The line starts from the catenary profile
// (straight without gravity or for a floating line) and is put into the static equilibrium of its
// discretization when the run starts (static_start, equilibrium()). A line state that is not finite after a step (an
// unstable explicit scheme) stops the run with the name of the line. Channels: end forces on the
// attached bodies, the largest segment tension and the tensions of the end segments (Ta, Tb), with BA.
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
#include"coral/mooring/coral_seabed.h"
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
        int point = -1;              // >= 0: a free point of the mooring system (body and p unused)
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
        double mu = 0.0, v_eps = 0.01;     // seabed friction of the dynamic line across it, regularization velocity
        double mu_axial = -1.0;            // seabed friction along the line (< 0: mu)
        double rho = 1025.0;
        double EI = 0.0;             // bending stiffness of the dynamic line [N m^2] (0: none, a cable or chain)
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

    // junctions of a line of sections (in order from a to b) between pa and pb: on the catenary of the whole line
    // (mean weight, series stiffness), the start of the junction points of a multi-section line
    std::vector<Vec3> mooring_junctions(const Vec3 &pa, const Vec3 &pb, const std::vector<mooring_props> &sections,
                                        double z_seabed, double g);

    // a free point of a mooring system (MoorDyn's free points): a bridle node, a clump weight, a buoy, the
    // junction of two lines. Its mass, buoyancy (volume, by the submerged fraction of the volume-equivalent
    // sphere), drag 1/2 rho CdA |u - v| (u - v), added mass rho V Ca and Froude-Krylov force rho V (1 + Ca) du/dt,
    // seabed contact (penalty kb, damping cb on the contact area pi r^2, friction mu) act together with the
    // ends of the attached lines (their tension, the weight, drag and added mass of the end half-segments).
    struct mooring_point_props
    {
        double mass = 0.0;           // [kg]
        double volume = 0.0;         // displaced volume [m^3]
        double CdA = 0.0;            // drag coefficient times area [m^2]
        double Ca = 0.0;             // added mass coefficient (of rho V)
        double radius = -1.0;        // contact radius on the seabed (< 0: of the volume-equivalent sphere, at least 0.05 m)
        double kb = 3.0e6, cb = 3.0e5, mu = 0.0, v_eps = 0.01;
    };

    // failure of a line: at a time or when the tension of an end segment (with BA) reaches T; the end lets go
    // of its body or point and becomes a free end (a point without mass of its own)
    struct mooring_break
    {
        double t = std::numeric_limits<double>::infinity();
        double T = std::numeric_limits<double>::infinity();
        int end = -1;                // 0: a, 1: b, -1: b for the time, the end that reaches T for the tension
    };

    // REEF3D::CORAL: a mooring system: dynamic lumped-mass lines (MoorDyn type, see above) whose ends sit on
    // bodies (or the ground) or on free points of the system. The inner nodes of all lines and the free
    // points are advanced together with explicit RK4 sub-steps inside the step of the multibody system;
    // the ends on bodies move on Hermite curves between the converged and the trial states of their
    // attachment points; the bodies carry the end nodes (their mass and added mass times the attachment's
    // acceleration). A system of one line without points is the mooring_line below.
    class mooring_system : public load
    {
    public:
        mooring_system(const model&, std::shared_ptr<const wave_kinematics> waves=nullptr, const Vec3 &current=Vec3::Zero(),
                       const std::string &name="mooring");

        // a free point at the global position p over a flat seabed at z_seabed (-1e30: none); returns its index
        // for mooring_end::point
        int add_point(const std::string &name, const Vec3 &p, const mooring_point_props&, double z_seabed=-1.0e30);
        // a line between a and b (an end on a body or the ground: body, global position; on a point:
        // point); z_seabed: height of a flat seabed under it (-1e30: none); returns the line index
        int add_line(const model&, const std::string &name, const mooring_end &a, const mooring_end &b, const mooring_props&,
                     int nseg, double z_seabed=-1.0e30);
        // sub-step of the explicit scheme (0: automatic, the smallest stable step of the lines and points)
        void set_dt_line(double dt);
        // the seabed of the lines and points that have one (a flat height above -1e30) is the surface (plane,
        // grid) instead of their flat heights; fluid_bed: the bed of the resolved fluid at the nodes where the
        // wave kinematics know it (wave_kinematics::bed_at), else the surface or the flat height (the lines and
        // points there and those added later); set_line_flat(k): line k keeps its flat seabed
        void set_seabed(std::shared_ptr<const seabed_surface> s, bool fluid_bed=false);
        void set_line_flat(int k) {lines_[k].surf = lines_[k].fbed = false;}
        // bending (EI > 0): end e of line k clamped to its body (or the ground) along the global direction d
        // into the line (zero: the line's start tangent there); the other ends are pinned
        void set_clamp(int k, int e, const model&, const Vec3 &d=Vec3::Zero());
        void set_point_flat(int j) {points_[j].surf = points_[j].fbed = false;}
        // failure of line k (checked at each converged step; channel <line>_broken)
        void set_break(int k, const mooring_break &b) {lines_[k].brk = b; lines_[k].has_break = true;}
        bool line_broken(int k) const {return lines_[k].broken >= 0;}
        double line_break_time(int k) const {return lines_[k].t_broken;}

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override;
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        void kinematics_points(const model&, const void *source, std::vector<Vec3>&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;
        // restart from a checkpoint (coral_state.h)
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;

        int nlines() const {return (int)lines_.size();}
        int npoints() const {return (int)points_.size();}
        const std::string &line_name(int k) const {return lines_[k].name;}
        const std::string &point_name(int j) const {return points_[j].name;}
        int find_point(const std::string&) const;                // -1 if none
        int line_nseg(int k) const {return lines_[k].N;}
        // converged node positions of line k (0 .. N)
        std::vector<Vec3> line_nodes(int k) const;
        Vec3 point_position(int j) const {return r_[points_[j].idx];}
        Vec3 point_velocity(int j) const {return v_[points_[j].idx];}
        double line_tension(int k, int seg) const;               // converged segment tension (with BA)
        // forces of the converged line k on its ends (on bodies with the inertia of the end nodes)
        void line_end_forces(int k, Vec3 &Fa, Vec3 &Fb) const {Fa = lines_[k].Fc[0]; Fb = lines_[k].Fc[1];}
        double dt_line() const {return dt_;}

        // relax the system with the body ends fixed (static shape), extra damping, before the run, in the
        // water at rest
        void settle(const model&, double t_settle);
        // wave kinematics at every Runge-Kutta stage (true) or once per step at the converged nodes (default)
        void set_exact_kinematics(bool b) {exact_kin_ = b;}
        // static equilibrium of the discretized lines and the free points between the current body ends in
        // the water at rest at t (damped Newton on the inner nodes and the points, Jacobian by coloured finite
        // differences: block-tridiagonal per line without points, sparse with points); true if converged (else
        // the shape is kept). initialize() calls it once (static_start, default), so the system starts at rest
        // instead of from the continuous catenaries, whose chords are shorter than the segments (a slack line
        // that drops and snaps)
        bool equilibrium(const model&, double t);
        void set_static_start(bool b) {static_start_ = b;}

    protected:
        struct state {std::vector<Vec3> r, v;};
        struct line_t
        {
            std::string name;
            mooring_props p;
            int N = 0;
            double l = 0.0;            // unstretched segment length
            double zbed = -1.0e30;     // flat seabed
            bool surf = false, fbed = false;               // seabed: the surface, the fluid's bed
            int off = 0;               // first node in the system's node arrays (N + 1 nodes)
            int pt[2] = {-1, -1};      // end on a free point (index), else -1
            int att[2] = {-1, -1};     // end on a body or the ground: attachment index
            double dt = 0.0;
            mooring_break brk;
            bool has_break = false;
            int broken = -1;           // the end that let go, -1: intact
            double t_broken = std::numeric_limits<double>::quiet_NaN();
            std::vector<Vec3> Tc;      // converged segment tensions
            Vec3 Fc[2] = {Vec3::Zero(), Vec3::Zero()};      // converged forces on the ends
            mutable std::vector<Vec3> T;                    // work: segment tensions
            mutable Vec3 F[2];                              // work: end forces of the last accel
            mutable Mat3 M[2];                              // work: end masses of the last accel
            bool clamp[2] = {false, false};                 // bending: end clamped to its body
            mutable Vec3 cd[2];                             // work: clamp directions (global) of the stage
            mutable Vec3 Mb[2] = {Vec3::Zero(), Vec3::Zero()};   // work: bending moments on the bodies
            mutable std::vector<Vec3> B;                    // work: bending forces
        };
        struct point_t
        {
            std::string name;
            mooring_point_props p;
            int idx = 0;               // node of the point in the node arrays
            double zbed = -1.0e30;     // flat seabed
            bool surf = false, fbed = false;
            std::vector<std::pair<int,int>> ends;           // attached (line, end)
            bool free_end = false;     // the free end of a broken line (no channels)
        };
        struct attach_t
        {
            int body = -1;
            Vec3 s = Vec3::Zero();     // in the body frame
            int line = 0, end = 0;
            int idx = 0;               // node in the node arrays
            Vec3 r = Vec3::Zero(), v = Vec3::Zero();        // converged
            Vec3 F = Vec3::Zero();     // converged force on the body (with the end inertia)
            Vec3 ds = Vec3::Zero();    // clamp direction in the body frame (zero: pinned)
            Vec3 d = Vec3::Zero();     // converged clamp direction (global)
        };
        struct end_motion {Vec3 r0, v0, r1, v1; Vec3 d0 = Vec3::Zero(), d1 = Vec3::Zero();};   // Hermite data of a body end over the step, clamp directions
        // wave kinematics at the converged nodes at time t (gen: the converged state they belong to)
        struct node_kin {double t = 0.0; long gen = -1; std::vector<double> eta, bed; std::vector<Vec3> u, du;};

        const node_kin &kinematics(int slot, double t) const;
        void water(int i, const Vec3 &x, double t, double &eta, Vec3 &u, Vec3 &du) const;
        // the seabed under node i at x (surface or fluid bed): height and unit normal; false: the flat seabed
        bool seabed_at(int i, const Vec3 &x, bool surf, bool fbed, double &zb, Vec3 &n) const;
        // accelerations of the inner nodes of line k for the state (its end nodes set), end forces and masses
        void line_accel(const line_t&, const state&, double t, std::vector<Vec3> &a, double extra) const;
        // accelerations of all inner nodes and points; the ends of the lines in s must be set
        void accel(const state&, double t, std::vector<Vec3> &a, double extra) const;
        // the line end nodes on points take the points' state
        void point_ends(state&) const;
        // from the converged state at t0 to t1 with the body ends on Hermite curves (motion per attachment)
        void integrate(state&, const std::vector<end_motion>&, double t0, double t1, double extra, bool at_rest=false) const;
        // keep the segment tensions and end forces of the converged state; throws if it is not finite
        void keep_converged();
        void update_dt();
        bool line_catenary(const line_t&, const Vec3 &pa, const Vec3 &pb, double g, Vec3 &Fa, Vec3 &Fb, std::vector<Vec3> *nodes) const;
        // quasi-static positions of the points on the catenaries of their lines (start of equilibrium_system)
        bool point_statics(const model&);
        // node layout: the lines' nodes (given), then the points (given)
        void layout(const std::vector<Vec3> &line_r, const std::vector<Vec3> &point_r);
        // end e of line k lets go: a new free point at the end's converged position and velocity
        void break_line(int k, int e, double t);
        bool equilibrium_lines(const model&, double t);     // each line on its own (no points)
        double roundoff(const line_t&) const;              // round-off of the accelerations of a line's nodes
        bool equilibrium_system(const model&, double t);    // lines and points together

        std::string name_;
        std::shared_ptr<const wave_kinematics> waves_;
        Vec3 current_;
        Vec3 g_;
        double dt_ = 0.0, dt_set_ = 0.0;
        std::shared_ptr<const seabed_surface> sb_;
        bool fbed_ = false;
        std::vector<line_t> lines_;
        std::vector<point_t> points_;
        std::vector<attach_t> att_;
        std::vector<int> free_;        // the nodes advanced by the scheme: inner nodes of the lines, then the points

        // converged state: all nodes at t_
        std::vector<Vec3> r_, v_;
        double t_ = 0.0;
        long gen_ = 0;                 // changes with the converged state
        bool static_start_ = true;     // initialize(): static equilibrium first (once)
        bool started_ = false;

        // wave kinematics: per step at the converged nodes (slots: start, end of the step) or exact
        bool exact_kin_ = false;
        mutable node_kin kin_[2];
        mutable int kmode_ = 0;        // 0: no waves, 1: exact, 2: linear between kin_[0] and kin_[1], 3: at rest (kin_[0] surface)
        mutable double kt0_ = 0.0, kt1_ = 0.0;
        mutable double creg_ = 0.0;    // static solve: compression stiffness of slack segments (fraction of EA)

        // work arrays
        mutable std::vector<Vec3> w_[4];
        mutable state ws_;

        // trial cache: time, positions and velocities of the attachments; forces and end masses
        mutable double ct_ = -1.0e300;
        mutable std::vector<Vec3> cr_, cv_, cF_, cMb_;
        mutable std::vector<Mat3> cM_;
        mutable std::vector<Vec3> cd_;     // clamp directions of the cached trial
        bool clamped_ = false;             // some end is clamped (bending)
        mutable bool cvalid_ = false;
    };

    // a single dynamic line between two attachment points: a mooring system of one line (the interface of
    // the line before the systems)
    class mooring_line : public mooring_system
    {
    public:
        mooring_line(const model&, const mooring_end &a, const mooring_end &b, const mooring_props&, int nseg,
                     std::shared_ptr<const wave_kinematics> waves=nullptr, const Vec3 &current=Vec3::Zero(),
                     double z_seabed=-1.0e30, double dt_line=0.0, const std::string &name="line");

        int nseg() const {return lines_[0].N;}
        std::vector<Vec3> nodes() const {return line_nodes(0);}  // converged node positions (0 .. N)
        double tension(int seg) const {return line_tension(0, seg);}
        // end forces of the converged line on the attached bodies (with the inertia of the end nodes)
        void end_forces(Vec3 &Fa, Vec3 &Fb) const {line_end_forces(0, Fa, Fb);}
    };
}

#endif
