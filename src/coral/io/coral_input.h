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
//   jforce <name> joint <joint> file <t F table> [scale <s>] [ramp <t>] [periodic] [rotational]
//        (prescribed force along the joint axis, or torque about it, interpolated in the table, 0 outside
//        it unless periodic; system identification with a multisine or chirp force, coral/pto/coral_joint_force.h)
//   spring <name> a <body> <x y z> b <body> <x y z> k <k> [c <c>] [L0 <L>] [T0 <T>] [tension_only]
//   linear6 <name> body <b> point <x y z> [k <k1..k6>] [c <c1..c6>]       (diagonal stiffness, damping)
//        [K <36 values>] [C <36 values>]           (full 6x6 matrices, row by row, e.g. coupled springs)
//   force <name> body <b> [point <x y z>] [F <fx fy fz>] [M <mx my mz>] [ramp <t>] [until <t>]
//        (constant load; until: released at that time, e.g. a free decay from the equilibrium under the load with static)
//   drag <name> body <b> [point <x y z>] (Q <q1..q6> | Cd <Cd> area <A> [axis <ex ey ez>]) [rho <rho>]
//        (quadratic damping F_i = -Q_i |u_i| u_i of the point velocity (default point: the CoG) and the
//        angular velocity, global axes; with Cd and area along the axis (default z) Q = rho Cd A / 2,
//        e.g. the viscous drag of a heaving float that a potential-flow fluid (FNPF) leaves out;
//        channels <name>_Fx .. _Mz, _P; coral/hydro/coral_body_drag.h)
//   line_type <name> d <d> mass <m/L> (EA <EA> | EA_curve <file>) [BA <BA>] [CB <CB>] [Cdn <c>] [Cdt <c>] [Can <c>] [Cat <c>]
//        [kb <k>] [cb <c>] [mu <mu>] [mu_axial <mu>] [mu_static <r>] [ks <k>] [cs <c>] [v_eps <v>] [EI <EI>] [BA_slack]
//        (properties of mooring lines for 'type' and 'section'; EA_curve: nonlinear elasticity of the dynamic line,
//        a file of 'strain tension' per line (strains increasing from >= 0, '#' comments), linear in between and
//        beyond the last point with its slope (MoorDyn's stress-strain table, which holds the last tension there);
//        the quasi-static line and the catenary start take the secant stiffness of the last point; BA < 0: -BA is
//        the damping ratio of the axial mode of a segment, BA = -BA l sqrt(EA m/L), l = L/nseg (MoorDyn's BA/-zeta))
//   mooring <name> <quasi_static|dynamic> a <body> <x y z> | a <point>  b <body> <x y z> | b <point>
//        (L <L> d <d> mass <m/L> EA <EA> | type <line type> L <L> | section <line type> <L> <nseg> [section ...])
//        [BA <BA>] [CB <CB>] [Cdn <c>] [Cdt <c>] [Can <c>] [Cat <c>] [kb <k>] [cb <c>] [mu <mu>]
//        [mu_axial <mu>] [mu_static <r>] [ks <k>] [cs <c>] [v_eps <v>] [BA_slack] (over those of the type)
//        [nseg <n>] [seabed <z> | no_seabed] [dt_line <dt>] [settle <T>] [current <ux uy uz>] [exact_kinematics]
//        [no_static] [break_time <t>] [break_tension <T>] [break_end <a|b>] [EI <EI>] [clamp_a [<dx dy dz>]] [clamp_b [...]]
//        [EA_curve <file>] [current_profile <file>] [outputs <p|t|pt>]
//        (current_profile: a file of 'z ux uy uz' per line, z increasing, linear in z, the end values outside, added to
//        'current' (MoorDyn's steady current profile); the lines of a mooring system take that of its first line;
//        outputs: channels of the dynamic line at its nodes, p the node positions <name>_n<i>_x .. _z (i = 0 .. nseg,
//        from end a), t the segment tensions with the axial damping <name>_T<j> (j = 0 .. nseg-1))
//        (sections: a line of several properties in order from a to b, e.g. chain-wire-chain: lines <name>.1 ..
//        <name>.n joined at points <name>.p1 .. (no mass or volume of their own), started on the catenary of the
//        whole line; dynamic lines only; failure: the line lets go of its end at break_time or when the tension of an
//        end segment reaches break_tension (break_end: that end only; default: b at the time, the end that reaches
//        the tension), the end then moves free; channel <name>_broken; EI: bending stiffness of a cable or riser
//        (dynamic lines), the ends pinned unless clamped to their body or the ground along a direction into the
//        line (default: the start tangent), the body then takes the bending moment)
//        (seabed: default the z of the ground end, the lower one if both ends are on the ground; a line between
//        two bodies (shared line, bridle) has no seabed unless given; the dynamic line takes the wave kinematics at
//        its converged nodes at the start and end of a step, linear in between (exact_kinematics: at every
//        Runge-Kutta stage, costly with many wave components); the dynamic line starts in the static
//        equilibrium of its segments between the ends at the start (no_static: from the continuous catenary,
//        whose chords are shorter than the segments); settle: relaxed in the water at rest before; seabed
//        friction of the dynamic line (default none): kinetic coefficients mu across and mu_axial (default mu) along
//        the line (elliptic limit), static ones mu_static times them; each node in contact sticks to an anchor
//        through a shear spring ks [Pa/m] and damper cs [Pa s/m] per contact area (default kb, cb) until the static
//        limit, then slides at the kinetic one, so a line stays where it lies under a load below the static limit
//        (on a slope across it, under a steady pull); the static start and settling are frictionless, the nodes stick
//        where they lie when the run starts. v_eps > 0: the regularized law of 0034 instead, mu Fn tanh(v / v_eps)
//        without static friction (a line creeps under any steady load). For example studless chain across 1.0,
//        along 0.69, mu_static 1.33 (the VolturnUS-S baseline of NREL's MoorDyn seabed friction study, 2023); use
//        values for the soil where known;
//        BA_slack: the axial damping BA in every segment of the dynamic line, also slack and compressed ones,
//        T = max(0, EA e) + BA de/dt (MoorDyn's and MoorDyn-F's law: a slack line resists folding up, snap loads
//        and a falling line as in MoorDyn); default: max(0, EA e + BA de/dt) while a segment is taut, nothing
//        when slack (a chain or rope never pushes);
//        channels <name>_Fa_x .. _Fb_z (forces on the bodies), _Tmax, _Ta, _Tb (end segment tensions, negative
//        for a segment that pushes with BA_slack); the quasi-static line has a seabed portion
//        with friction CB only at a ground end on the seabed, else it hangs free between its ends; the names of lines and other loads must not be keywords of the line, e.g. not "d" or "L"; the
//        dynamic line takes the free surface and kinematics of the waves of coral.txt, of the resolved fluid
//        with fluid_kinematics, or else still water at z = 0 of CORAL's frame: with REEF3D set
//        fluid_coupling origin 0 0 <F 60> or fluid_kinematics, CORAL warns otherwise)
//   bathymetry plane <z0> [<sx> <sy>] | grid <file> | fluid
//        (the seabed of the dynamic mooring lines and points that have one, instead of their flat heights, except those
//        with 'seabed <z>': the plane z0 + sx x + sy y, a grid of heights (header 'nx ny x0 y0 dx dy', then nx*ny heights,
//        x fastest; bilinear, its edge outside), or the bed of the resolved fluid at the nodes (FNPF, NHFLOW, with
//        fluid_kinematics; elsewhere the flat height); contact along the normal of the plane or grid, vertical on the
//        fluid's bed; the quasi-static lines keep their flat seabed)
//   point <name> <x y z> [mass <m>] [volume <V>] [CdA <A>] [Ca <c>] [radius <r>] [kb <k>] [cb <c>] [mu <mu>] [mu_static <r>]
//        [ks <k>] [cs <c>] [v_eps <v>]
//        [seabed <z> | no_seabed]
//        (a free point of a mooring system: bridle node, clump weight, buoy, junction of two lines; dynamic lines
//        end on it with 'a <point>' or 'b <point>' (no coordinates); weight, buoyancy of the submerged volume,
//        drag 1/2 rho CdA |u - v| (u - v), added mass rho V Ca, seabed contact of a sphere of the radius (default:
//        volume-equivalent; without volume or radius only the line ends touch the seabed; friction mu, mu_static, ks, cs
//        as for the lines, isotropic); the lines connected through points form one mooring system (advanced together),
//        whose seabed defaults to the lowest ground end of its lines; channels <point>_x, _y, _z)
//   morison <name> body <b> p1 <x y z> p2 <x y z> D <D> [Cd <c>] [Ca <c>] [Cm <c>] [n <n>] [current <u v w>]
//   strip_member <name> body <b> p1 <x y z> p2 <x y z> D <D> [Cd <c>] [Ca <c>] [Cm <c>] [n <n>] [end1] [end2]
//        [Ca_end <c>] [Cd_end <c>] [no_buoyancy] [current <u v w>] [exact_kinematics]
//        (strip theory of a rigid floater with nonlinear hydrostatics: buoyancy of the wet part of each strip below the
//        instantaneous elevation, Morison normal to the member, and at closed ends below the surface (end1, end2) the
//        waves' dynamic pressure, added mass Ca_end rho D^3/6 and drag along the axis; channels <name>_Fx .. _Mz
//        (moment about the body's CoG), _V (submerged volume); the wave kinematics at the undisplaced strips once per
//        time unless exact_kinematics; coral/hydro/coral_strip_hydro.h)
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
//   rho_air <rho>                                air density for the aerodynamic loads (1.225)
//   wind steady U <U> [z_ref <z>] [direction <deg>] [shear <alpha> | log <z0>] [veer <deg/m>] [upflow <deg>]
//        [origin <x y>] [ramp <t>]
//   wind file <uniform wind file> [z_ref <z>] [L_ref <L>] [shear ...] [veer ...] [origin <x y>] [ramp <t>]
//        (U at the height z_ref (default: z_hub of wind iec), direction of propagation counterclockwise from +x, power-law
//        or logarithmic shear, veer: turning of the direction with height, origin: the hub's horizontal position (lateral
//        shear, events, turbulence); the file: InflowWind's uniform wind, columns t U direction (deg, positive clockwise
//        from above) [w] [horizontal shear] [shear exponent] [linear vertical shear] [gust], the linear shears per L_ref
//        (default the rotor diameter of wind iec), linear in time; coral/wind/coral_wind.h)
//   wind iec <I|II|III> <A+|A|B|C> z_hub <z> D <D>    IEC 61400-1 ed. 3 turbine class: events and turbulence
//   wind event <eog|edc+|edc-|ecd+|ecd-|ews_v+|ews_v-|ews_h+|ews_h-> t0 <t>
//        (extreme operating gust, direction change, coherent gust with direction change, vertical or horizontal shear,
//        from t0, at the hub speed of the steady wind; normal turbulence model sigma_1)
//   wind turbulence kaimal [ntm|etm|ewm] [ny <n>] [nz <n>] [width <w>] [height <h>] [T <T>] [dt <dt>] [seed <s>]
//        [coherence_vw] [scale_hub] [write <file.bts>]
//        (IEC Kaimal turbulence with the exponential coherence (u; v and w too with coherence_vw) on a grid around the
//        hub, default 15 x 15 points over 1.1 D, T end + 10 s, dt 0.05 s, rounded up to 2^n levels (periodic); frozen and
//        convected with the hub speed; write: the full field (mean + turbulence) as a TurbSim .bts for other codes;
//        coral/wind/coral_wind_grid.h)
//   wind turbsim <file.bts> [direction <deg>] [origin <x y>] [x0 <xi0>]
//        (TurbSim full field, frozen: a point at the downwind distance xi from the origin sees the grid at
//        t + (xi0 - xi)/u_hub; outside the grid its edge, below it the tower points)
//   wind_probe <name> <x y z> [body <b>]          (output channels <name>_u _v _w _U of the inflow at the point)
//   rotor <name> hub <x y z> axis <ax ay az> R_hub <r> blade <AeroDyn blade file> airfoils <n> <file 1> .. <file n>
//        [blades <B>] [structure <ElastoDyn blade file>] [hub_mass <m>] [hub_inertia <J> [<J_across>]] [precone <deg>]
//        [azimuth <deg>] [ccw] [nacelle <body>] [rpm <rpm>] [pitch <deg> | pitch_file <t deg table>]
//        [gearbox <N>] [gen_inertia <J_hss>] [shaft <k> <c>] [efficiency <eta>]
//        [generator none | constant <T> | slip <slope> <rpm_sync> | k_omega2 <K> [<T_max>]]
//        [no_skew] [no_tip_loss] [no_hub_loss] [no_tangential] [no_induction] [no_drag_induction]
//        [tower <x y z base> <x y z top> <D_base> <D_top> [tower_body <b>]]
//        [pitch_joints [pitch_freq <Hz>] [pitch_zeta <z>]] [blade_elements <n> [blade_beta <s>] [ea_per_mass <EA/m>] [gj_factor <f>]]
//        (wind turbine rotor: hub (body <name>_hub, <name>_gen with the shaft torsion) turning in the revolute joint
//        <name>_shaft on the nacelle (default the ground) about the axis
//        pointing downwind, clockwise looking downwind unless ccw; blade 1 up at azimuth 0; the stations of AeroDyn's
//        blade file at R_hub + BlSpn from the hub centre, airfoils in the order of its BlAFID; the blade mass from
//        ElastoDyn's blade file (BMassDen); precone positive with the tips downwind; generator torque on the high-speed
//        shaft (N m, rad/s; rpm_sync of the high-speed shaft), shaft torsion k [N m/rad], c [N m s/rad] on the
//        low-speed shaft with the generator on a body of its own; quasi-steady blade element momentum with
//        Prandtl's tip and hub loss, Buhl's correction and the skewed wake; tower potential flow;
//        channels <name>_azimuth _rpm _T _Q _P _Cp _Ct _lambda _U _yaw _pitch _b<k>_Moop _b<k>_Mip and
//        <name>_dt_rpm_rotor _rpm_gen _T_gen _P_gen _T_shaft; blades rigid on the hub unless: pitch_joints: each blade
//        on a pitch bearing (body <name>_pitch<k>, revolute joint <name>_pitch<k> on the hub, positive towards feather)
//        driven to the pitch command by a second-order actuator (default 2 Hz, damping 0.7; channels <name>_pitch_cmd,
//        _pitch_b<k>, _pitch_b<k>_T); blade_elements n: flexible blades of n co-rotational beam elements from the
//        structural table (structure is then required; EA = ea_per_mass m (2e7), GJ = gj_factor EI_flap (0.3), damping
//        blade_beta (0.0022 s)), nodes <name>_b<k>e<i>:1; coral/aero/coral_rotor.h)
//   controller <name> rotor <rotor> [tune P_rated <W> rpm_rated <rpm> [omega_n <rad/s>] [zeta <z>] [write <file>]]
//        [K <N m s^2>] [rpm_cutin <rpm>] [rpm_rgn2 <rpm>] [rpm_rated_torque <rpm>] [P_rated <W>] [slip <%>] [T_max <N m>]
//        [T_rate <N m/s>] [pitch_rgn3 <deg>] [rpm_ref <rpm>] [KP <s>] [KI <->] [theta_K <deg>] [theta_min <deg>]
//        [theta_max <deg>] [theta_rate <deg/s>] [corner <Hz>] [dt <s>] [no_pitch] [no_torque]
//        (wind turbine controller of the NREL 5 MW baseline type on the filtered generator speed: variable-speed
//        torque (regions 1, 1.5, 2 k omega^2, 2.5, 3 constant power) and gain-scheduled PI collective pitch; speeds and
//        torque on the generator (high-speed) side, K in N m s^2/rad^2; the defaults are the NREL 5 MW's (DISCON); it sets
//        the rotor's generator torque and the pitch command of its pitch actuators (or the pitch of rigid blades without
//        pitch bearings) once per step; tune: the gains for this rotor from CORAL's steady BEM at the rated power and the
//        rated rotor speed (K from Cp_max and lambda_opt, the PI gains from dP/dtheta for a response at omega_n
//        (0.6 rad/s) and zeta (0.7) with the rotor's and generator's inertia, theta_K from the fit of dP/dtheta,
//        rpm_ref = N rpm_rated (the rotor speed), rpm_rated_torque (region 3 of the torque law) 0.99 of it, rpm_cutin, rpm_rgn2 0.57 and 0.74 of rpm_ref
//        and T_max 1.1 rated torque as the NREL 5 MW; K Omega_rated^3 must stay below P_rated: the baseline torque law
//        reaches the rated power at the rated speed, not before), given values override; write: the steady operating curve and the
//        tuning; channels <name>_wf (rpm), _T, _pitch (deg), _region; coral/control/coral_turbine_controller.h)
//   static                                       static equilibrium before the time integration
//   modes <n> [damped]                           linear analysis: the first n modes (modes.dat, mode_shapes.dat)
//   output vtp <dt>                              VTP of the structure (coral_vtp/)

#ifndef CORAL_INPUT_H_
#define CORAL_INPUT_H_

#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/hydro/coral_waves.h"
#include"coral/wind/coral_wind.h"
#include"coral/aero/coral_rotor.h"
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
        std::shared_ptr<wind_inflow> wind;
        std::map<std::string, rotor_handles> rotors;       // rotor <name>
        double rho_air = 1.225;
        std::map<std::string, joint*> joints;
        double rho = 1025.0;

        // bodies driven by the resolved fluid of REEF3D (fluid_body) and their load
        struct fluid_body_spec {int body = -1; int n = 0; bool own_mass = false; bool one_way = false; bool excitation = false;
                               Vec6 added_mass = Vec6::Zero();};
        std::vector<fluid_body_spec> fluid_bodies;
        fluid_load *fluid = nullptr;
        std::shared_ptr<sampled_kinematics> fluid_kin;     // fluid_kinematics
        bool fluid_bed = false;                             // bathymetry fluid: the mooring lines on the fluid's bed
        struct fluid_member_spec {member_fluid_load *load = nullptr; int n = 0;};
        std::vector<fluid_member_spec> fluid_members;       // fluid_member
        Vec3 fluid_origin = Vec3::Zero();
        std::string fluid_record;                           // fluid_coupling record <file>
        int fluid_order = 2;
        bool fluid_predictor = true;

    private:
        std::map<std::string, rotor_spec> rotor_specs_;
        void parse(const std::vector<std::pair<int, std::vector<std::string>>> &lines, const std::string &dir);
    };
}

#endif
