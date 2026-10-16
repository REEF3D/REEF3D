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
//        [latching <t_hold> <c_latch>] [declutching <t_hold>]
//   spring <name> a <body> <x y z> b <body> <x y z> k <k> [c <c>] [L0 <L>] [T0 <T>] [tension_only]
//   linear6 <name> body <b> point <x y z> [k <k1..k6>] [c <c1..c6>]       (diagonal stiffness, damping)
//   mooring <name> <quasi_static|dynamic> a <body> <x y z> b <body> <x y z> L <L> d <d> mass <m/L>
//        EA <EA> [BA <BA>] [CB <CB>] [Cdn <c>] [Cdt <c>] [Can <c>] [Cat <c>] [kb <k>] [cb <c>] [mu <mu>]
//        [nseg <n>] [seabed <z>] [dt_line <dt>] [settle <T>] [current <ux uy uz>]
//   morison <name> body <b> p1 <x y z> p2 <x y z> D <D> [Cd <c>] [Ca <c>] [Cm <c>] [n <n>] [current <u v w>]
//   section <name> tube D <D> t <t> E <E> G <G> rho <rho>
//   section <name> EA <EA> EIy <EI> EIz <EI> GJ <GJ> m <m/L> [GAy <GA>] [GAz <GA>] [rIx <r>] [rIy <r>] [rIz <r>] [D <D>]
//   member <name> p1 <x y z> p2 <x y z> n <n> section <s> [section2 <s>] [orient <x y z>] [beta <b>]
//        [end1 <body|ground>] [end2 <body|ground>]
//        (co-rotational beam elements on nodes "<name>:0" .. "<name>:n", usable as bodies by the items
//        below; members meeting at a point share the node; end1/end2 put an end on a body or the ground)
//   point_mass <body> mass <m> [inertia <Ixx Iyy Izz>]       (added to a body, e.g. a node)
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

    private:
        void parse(const std::vector<std::pair<int, std::vector<std::string>>> &lines, const std::string &dir);
    };
}

#endif
