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
// REEF3D::CORAL: the CORAL input file (model and run), read into a coral_case.
//
// One item per line, '#' starts a comment, a line ending in '\' continues on the next line.
// Keywords followed by key-value pairs; bodies and joints are referred to by name, "ground" is the
// fixed frame. Angles in degrees, SI units, z up, still water level z = 0.
//
//   gravity <gx> <gy> <gz>                        default 0 0 -9.81 (gravity on all bodies)
//   rho <rho>                                     water density for Morison and mooring (1025)
//   time dt <h> end <t> [start <t0>] [rho_inf <r>] [tol <tol>] [full_newton]
//   output dt <h> [file <name>]
//   body <name> mass <m> [inertia <Ixx Iyy Izz> [<Ixy Ixz Iyz>]] [cog <x y z>] [rpy <r p y>]
//        [velocity <vx vy vz>] [omega <wx wy wz>] [offset <dx dy dz>]
//        (inertia about the CoG, body axes; the joints and loads are set up in the configuration given by
//        cog and rpy, the equilibrium; offset displaces the body afterwards, e.g. for a decay test)
//   joint <name> <fixed|spherical|revolute|prismatic|cylindrical|planar> <body a> <body b>
//        point <x y z> [axis <ax ay az>]
//   waves regular H <H> T <T> [direction <deg>] [phase <deg>]
//   waves component a <a> w <w> [direction <deg>] [phase <deg>]
//   waves <jonswap|pm> Hs <Hs> Tp <Tp> [gamma <g>] [wmin <w>] [wmax <w>] [n <n>] [seed <s>] [direction <deg>]
//   waves depth <d> | ramp <t>                    (apply to all components; depth <= 0: infinite)
//   hydro <name> <file <db> | wamit <root> ulen <L>> bodies <b1> [<b2> ...] [state_space]
//        [memory <T> <dt>] [no_radiation] [no_excitation] [no_hydrostatics] [no_buoyancy]
//        [damping <mode> <value>] [drag <mode> <value>] [volume <i> <V>] [ref <i> <x y z>]
//        (modes 1..6N of the database; damping linear, drag quadratic, diagonal entries)
//   pto <name> joint <joint> [rotational] [k <k>] [c <c>] [q0 <q>] [F0 <F>] [coulomb <Fc>]
//        [v_eps <v>] [F_max <F>] [P_max <P>] [stops <qmin> <qmax> <k_stop> [<c_stop>]]
//        [latching <t_hold> <c_latch> [<k_latch>]] [declutching <t_hold>]
//        (elements and controllers of the PTO library, coral/pto/coral_pto_element.h)
//   drive <name> joint <joint> sine <amplitude> <period> [phase <deg>] [ramp <t>] [q0 <q>] [rotational]
//        [k <k>] [c <c>]          (prescribed joint motion by a stiff spring-damper, default k 1e7, c 1e4;
//        forced oscillation in a resolved fluid, coral/pto/coral_joint_drive.h)
//   spring <name> a <body> <x y z> b <body> <x y z> k <k> [c <c>] [L0 <L>] [T0 <T>] [tension_only]
//   linear6 <name> body <b> point <x y z> [k <k1..k6>] [c <c1..c6>]       (diagonal stiffness, damping)
//        [K <36 values>] [C <36 values>]           (full 6x6 matrices, row by row, e.g. coupled springs)
//   force <name> body <b> [point <x y z>] [F <fx fy fz>] [M <mx my mz>] [ramp <t>]   (constant load)
//   drag <name> body <b> [point <x y z>] (Q <q1..q6> | Cd <Cd> area <A> [axis <ex ey ez>]) [rho <rho>]
//        (quadratic damping F_i = -Q_i |u_i| u_i of the point velocity (default point: the CoG) and the
//        angular velocity, global axes; with Cd and area along the axis (default z) Q = rho Cd A / 2,
//        e.g. the viscous drag of a heaving float that a potential-flow fluid (FNPF) leaves out;
//        channels <name>_Fx .. _Mz, _P; coral/hydro/coral_body_drag.h)
//   line_type <name> d <d> mass <m/L> EA <EA> [BA <BA>] [CB <CB>] [Cdn <c>] [Cdt <c>] [Can <c>] [Cat <c>] [kb <k>] [cb <c>]
//        [mu <mu>] [mu_axial <mu>] [v_eps <v>]       (properties of mooring lines for 'type' and 'section')
//   mooring <name> <quasi_static|dynamic> a <body> <x y z> | a <point>  b <body> <x y z> | b <point>
//        (L <L> d <d> mass <m/L> EA <EA> | type <line type> L <L> | section <line type> <L> <nseg> [section ...])
//        [BA <BA>] [CB <CB>] [Cdn <c>] [Cdt <c>] [Can <c>] [Cat <c>] [kb <k>] [cb <c>] [mu <mu>]
//        [mu_axial <mu>] [v_eps <v>] (over those of the type)
//        [nseg <n>] [seabed <z> | no_seabed] [dt_line <dt>] [settle <T>] [current <ux uy uz>] [exact_kinematics]
//        [no_static]
//        (sections: a line of several properties in order from a to b, e.g. chain-wire-chain: lines <name>.1 ..
//        <name>.n joined at points <name>.p1 .. (no mass or volume of their own), started on the catenary of the
//        whole line; dynamic lines only)
//        (seabed: default the z of the ground end, the lower one if both ends are on the ground; a line between
//        two bodies (shared line, bridle) has no seabed unless given; the dynamic line takes the wave kinematics at
//        its converged nodes at the start and end of a step, linear in between (exact_kinematics: at every
//        Runge-Kutta stage, costly with many wave components); the dynamic line starts in the static
//        equilibrium of its segments between the ends at the start (no_static: from the continuous catenary,
//        whose chords are shorter than the segments); settle: relaxed in the water at rest before; seabed
//        friction mu across and mu_axial (default mu) along the line, regularized below v_eps (0.01 m/s);
//        channels <name>_Fa_x .. _Fb_z (forces on the bodies), _Tmax, _Ta, _Tb (end segment tensions); the quasi-static line has a seabed portion
//        with friction CB only at a ground end on the seabed, else it hangs free between its ends; the names of lines and other loads must not be keywords of the line, e.g. not "d" or "L"; the
//        dynamic line takes the free surface and kinematics of the waves of coral.txt, of the resolved fluid
//        with fluid_kinematics, or else still water at z = 0 of CORAL's frame: with REEF3D set
//        fluid_coupling origin 0 0 <F 60> or fluid_kinematics, CORAL warns otherwise)
//   point <name> <x y z> [mass <m>] [volume <V>] [CdA <A>] [Ca <c>] [radius <r>] [kb <k>] [cb <c>] [mu <mu>] [v_eps <v>]
//        [seabed <z> | no_seabed]
//        (a free point of a mooring system: bridle node, clump weight, buoy, junction of two lines; dynamic lines
//        end on it with 'a <point>' or 'b <point>' (no coordinates); weight, buoyancy of the submerged volume,
//        drag 1/2 rho CdA |u - v| (u - v), added mass rho V Ca, seabed contact of a sphere of the radius (default:
//        volume-equivalent); the lines connected through points form one mooring system (advanced together),
//        whose seabed defaults to the lowest ground end of its lines; channels <point>_x, _y, _z)
//   morison <name> body <b> p1 <x y z> p2 <x y z> D <D> [Cd <c>] [Ca <c>] [Cm <c>] [n <n>] [current <u v w>]
//   morison <name> member <m> [D <D>] [Cd <c>] [Ca <c>] [Cm <c>] [n <strips per element>] [current <u v w>]
//        [no_buoyancy] [flooded <D_inner>]
//        (distributed on a frame or superelement member, buoyancy of the submerged strips; D defaults to
//        the section's; flooded: the inner water below z = 0 as node mass, frame members only)
//   section <name> tube D <D> t <t> E <E> G <G> rho <rho>
//   section <name> EA <EA> EIy <EI> EIz <EI> GJ <GJ> m <m/L> [GAy <GA>] [GAz <GA>] [rIx <r>] [rIy <r>] [rIz <r>] [D <D>]
//   member <name> p1 <x y z> p2 <x y z> n <n> section <s> [section2 <s>] [orient <x y z>] [beta <b>]
//        [end1 <body|ground>] [end2 <body|ground>]
//        (co-rotational beam elements on nodes "<name>:0" .. "<name>:n", usable as bodies by the items
//        below; members meeting at a point share the node; end1/end2 put an end on a body or the ground)
//   point_mass <body> mass <m> [inertia <Ixx Iyy Izz>] [at <x y z>]
//        (added to a body, e.g. a node; with at: rigidly attached at that global point)
//   superelement <name> members <m1> [<m2> ...] [modes <n>] [zeta <z>] [beta <b>]
//        (the listed members form a Craig-Bampton superelement instead of frame members: ends on the
//        ground clamped, ends on bodies are the interface; n fixed-interface modes (10), modal damping
//        ratio, stiffness-proportional damping of the interface part)
//   soil <name> mudline <z>
//   soil_layer <soil> sand top <d> bottom <d> phi <deg> k <N/m3> [gamma <N/m3>] [delta <deg>] [K <K>] [Nq <Nq>] [cyclic]
//   soil_layer <soil> clay top <d> bottom <d> su <Pa> [su_bottom <Pa>] [gamma <N/m3>] [eps50 <e>] [J <J>]
//        [residual <r>] [cyclic]
//   soil_layer <soil> linear top <d> bottom <d> kh <N/m2> [kv <N/m2>] [kq <N/m>]
//        (depths below the mudline; gamma: submerged unit weight; API p-y, t-z and Q-z curves, see
//        soil/coral_soil.h)
//   pile <frame member> soil <soil> D <D> [py] [tz] [qz] [linear]
//        (springs at the member's nodes below the mudline; all three kinds unless some are listed;
//        linear: initial slopes only)
//   fluid_body <body> [n <k>] [own_mass] [added_mass <a1 .. a6>] [one_way [excitation]]
//        (REEF3D with Z 40 1: the body is driven by CORAL, the resolved fluid acts on it; k: number of the
//        REEF3D 6DOF body (X 20 order, default: the order of the fluid_body lines); mass and inertia are
//        taken from the REEF3D body unless own_mass; FNPF, NHFLOW and CFD, see coral/reef3d.
//        added_mass: an artificial added mass (diagonal, surge sway heave roll pitch yaw, global axes about
//        the CoG) for a fluid without its own (NHFLOW, CFD): implicit in CORAL and compensated by the same mass
//        times the acceleration the fluid saw; it stabilises the staggered coupling of light bodies. Use it
//        only where needed: NHFLOW reacts to the body acceleration within a step with much less than the
//        added mass of the body's motion: half to once that added mass stabilised a box 9 times lighter
//        than it (validation 08) at a few percent error in the start transient, a quarter did not, and
//        twice made the 4 kg box unstable too; for bodies not much lighter than their added mass leave it
//        out (more accurate). CFD reacts within a step with about the full added mass: a 1.5 kg box in waves
//        (added mass about 4 kg) needed about once its added mass in heave and pitch (half of it was
//        unstable), a 0.3 kg box (added mass 0.37 kg) was stable from zero to once (validation 10).
//        one_way: REEF3D keeps its body fixed (X 11 0) or prescribed, only its loads act on the CORAL body,
//        at the REEF3D CoG, which may be any CORAL body (FNPF, NHFLOW and CFD); excitation: without the load of
//        the start (still water), e.g. with the hydrostatics and radiation of a hydro database)
//   fluid_member <name> member <m> [n <k>] [strips <ns>] [order <0|1|2>] [excitation]
//        (REEF3D with Z 40 1, FNPF, NHFLOW or CFD (X 60 1): the loads on the hull of the REEF3D 6DOF body k, kept fixed or
//        prescribed there, integrated in ns strips along the frame or superelement member m (default two per
//        element) and applied along it; order of the time polynomial (2); excitation: without the loads of
//        the start; see coral/hydro/coral_member_fluid_load.h)
//   fluid_kinematics [order <0|1|2>]
//        (REEF3D with Z 40 1, FNPF or NHFLOW: the Morison members (on bodies and frame/superelement members)
//        and the dynamic mooring lines take the velocity, acceleration and free surface of the resolved
//        fluid, sampled at their points once per fluid step and interpolated in time with this order
//        (default 2); the waves of coral.txt act outside the fluid domain and before the first sample)
//   fluid_coupling [order <0|1|2>] [no_predictor] [origin <x y z>] [record <file>]
//        (fluid loads in time over a CORAL step: polynomial of this order through the last loads of the
//        fluid, default 2; predictor (default): the load of the next fluid stage, computed at a predicted
//        pose, enters the CORAL step, which then interpolates; see coral/hydro/coral_fluid_load.h);
//        origin: the point of the REEF3D grid that is CORAL's origin (default 0 0 0; e.g. 0 0 <depth> for
//        CORAL's still water level z = 0 with REEF3D's bed at z = 0); record: a journal of all fluid input
//        (REEF3D_CORAL/<file>) for the standalone CORAL: CORAL coral.txt replay <file> (io/coral_coupling.h)
//   static                                       static equilibrium before the time integration
//   modes <n> [damped]                           linear analysis: the first n modes (modes.dat, mode_shapes.dat)
//   output vtp <dt>                              VTP of the structure (coral_vtp/)

#ifndef CORAL_INPUT_H_
#define CORAL_INPUT_H_

#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/hydro/coral_waves.h"
#include<map>
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    class joint;
    class fluid_load;
    class sampled_kinematics;
    class member_fluid_load;

    struct run_settings
    {
        double t_start = 0.0;
        double t_end = 10.0;
        double dt = 0.01;
        double out_dt = 0.0;          // 0: every step
        std::string out_file = "coral_results.dat";
        double vtp_dt = -1.0;         // < 0: no VTP output of the structure
        bool static_ini = false;      // static equilibrium before the time integration
        int modes = 0;                // number of modes of the linear analysis (0: none)
        bool modes_damped = false;
        timeint_settings ti;
    };

    class coral_case
    {
    public:
        // read the input file; throws std::runtime_error with the line number on errors
        void read(const std::string &file);
        // the same from a string (tests)
        void read_string(const std::string &text, const std::string &dir=".");

        model m;
        run_settings run;
        std::shared_ptr<wave_field> waves;
        std::map<std::string, joint*> joints;
        double rho = 1025.0;

        // bodies driven by the resolved fluid of REEF3D (fluid_body) and their load
        struct fluid_body_spec {int body = -1; int n = 0; bool own_mass = false; bool one_way = false; bool excitation = false;
                               Vec6 added_mass = Vec6::Zero();};
        std::vector<fluid_body_spec> fluid_bodies;
        fluid_load *fluid = nullptr;
        std::shared_ptr<sampled_kinematics> fluid_kin;     // fluid_kinematics
        struct fluid_member_spec {member_fluid_load *load = nullptr; int n = 0;};
        std::vector<fluid_member_spec> fluid_members;       // fluid_member
        Vec3 fluid_origin = Vec3::Zero();
        std::string fluid_record;                           // fluid_coupling record <file>
        int fluid_order = 2;
        bool fluid_predictor = true;

    private:
        void parse(const std::vector<std::pair<int, std::vector<std::string>>> &lines, const std::string &dir);
    };
}

#endif
