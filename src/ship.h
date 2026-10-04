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

#ifndef SHIP_H_
#define SHIP_H_

#include"6DOF_load.h"
#include"ship_models.h"
#include<vector>
#include<fstream>
#include<string>

class lexer;

using namespace std;

//  Ship module (X 350 1): semi-empirical hull loads, propulsion and steering of a 6DOF hull,
//  read from ship.dat in the case folder. Solver independent: the hull is the 6DOF surface
//  triangulation, the state that of the rigid-body core; the propeller acts on the fluid
//  through the actuator disk of the coupling (NHFLOW).
//
//  ship.dat (one keyword per line, # comments; all optional):
//
//  hull
//   lpp            L            length for the Reynolds number [m]     (default: waterline length)
//   wetted_surface S            [m^2]                                 (default: hull below the still water level)
//   form_factor    k            (1+k) C_F                              (default 0)
//   friction       0|1          ITTC-1957 frictional resistance        (default 1; 0 for CFD and with X 38 1)
//   viscosity      nu           kinematic viscosity [m^2/s]            (default W 2)
//   roll_damping   B44 B44q     linear [N m s] and quadratic [N m s^2]  (default 0 0)
//   crossflow      Cd           cross-flow drag coefficient            (default 0: off)
//   strips         n            strips for the cross-flow drag         (default 40)
//   thrust         T [x z]      constant thrust along the ship x-axis [N] at (x, z) relative
//                               to the CoG in the ship frame           (default 0)
//
//  propeller (body force, Hough-Ordway actuator disk)
//   propeller      x y z D      disk centre relative to the CoG [m], diameter [m]
//   propeller_hub  rh           hub radius / tip radius                (default 0.2)
//   propeller_kt   kt0 kt1 kt2  KT(J) = kt0 + kt1 J + kt2 J^2
//   propeller_kq   kq0 kq1 kq2  KQ(J) likewise
//   propeller_rps  n            revolutions per second                 (default 0)
//   propeller_sense s           +1 right-handed (clockwise seen from aft), -1 left (default 1)
//   propeller_thickness dx      axial extent of the actuator disk [m]  (default max(0.2 D, 4 dx))
//   propeller_inflow wake w     axial inflow Va = (1-w) u              (default, w = 0)
//   propeller_inflow sample d   Va from the fluid velocity on a ring d [m] ahead of the disk
//   propeller_source 0|1        actuator disk in the fluid; 0: thrust (1-t) T on the hull only
//                               (default 1 for NHFLOW, 0 otherwise)
//   thrust_deduction t          for propeller_source 0                 (default 0)
//
//  rudder (MMG standard method)
//   rudder         xR zR AR Lambda      position (xR < 0 aft) [m], area [m^2], aspect ratio
//   rudder_mmg     tR aH xH eps kappa lR gammaR    (defaults 0.39 0.3 -0.45 L 1.1 0.5 -0.9 L 0.4)
//   rudder_angle   fixed d              rudder angle [deg] (> 0: to starboard)
//   rudder_angle   zigzag d psi [t0]    zig-zag manoeuvre d/psi [deg], starting at t0 [s]
//   rudder_angle   autopilot psi Kp Kd [Ki]   heading autopilot, target heading [deg]
//   rudder_rate    r            maximum rudder rate [deg/s]            (default 15)
//   rudder_max     d            maximum rudder angle [deg]             (default 35)
//
//  The ship frame is the body frame of the 6DOF object (x forward, y to port, z up at the start);
//  headings psi are counter-clockwise from x.
//  Output: REEF3D_<model>_6DOF/REEF3D_ship_<n>.dat

class ship : public sixdof_load
{
public:
    
    ship(lexer*, int);
    virtual ~ship();
    
    void add_load(lexer*, const sixdof_rigidbody&, const sixdof_geometry&, sixdof_fluid*, double*) override;
    void actuator_disks(vector<sixdof_actuator_disk>&) const override;
    void print(lexer*) override;
    
private:
    
    void read(lexer*);
    void ini(lexer*, const sixdof_rigidbody&, const sixdof_geometry&);
    void steering(lexer*, double, double);
    double propeller_inflow(const sixdof_rigidbody&, sixdof_fluid*, const Eigen::Vector3d&, const Eigen::Vector3d&);
    
    const int id;
    bool initialized;
    
    // hull
    double lpp, S, k, nu, B44, B44q, Cd, thrust, xthrust, zthrust;
    int friction, nstrip;
    bool lpp_in, S_in;
    double xa, xf, zw;
    vector<double> xs, dx, T;
    
    // propeller
    bool prop;
    double xp, yp, zp, Dp, hub, kt[3], kq[3], nrps, thick, wake, sample_d, tded;
    int sense, inflow_mode, psource;
    double Va, J, KT, KQ, Tp, Qp;
    sixdof_actuator_disk disk;
    
    // rudder
    bool rud;
    ship_models::rudder_param rp;
    bool rmmg_in;
    int rmode;
    double rcmd, zz_d, zz_psi, zz_t0, ap_psi, ap_Kp, ap_Kd, ap_Ki;
    double rrate, rmax;
    double delta, tlast, psi_c, psi_last, psi0, eint;
    int zz_sign;
    bool zz_started;
    double alphaR, UR, FN, XR, YR, NR, KR;
    
    // loads of the last evaluation, ship frame
    double ub, vb, wb, pb, qb, rb_;
    double Re, CF, XF, Ycf, Ncf, Kroll;
    
    ofstream out;
};

#endif
