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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SHIP_H_
#define SHIP_H_

#include"6DOF_load.h"
#include<vector>
#include<fstream>
#include<string>

class lexer;

using namespace std;

//  Ship module (X 350 1): semi-empirical loads on a 6DOF hull that the hydrodynamic model
//  does not resolve, read from ship.dat in the case folder. Solver independent: the hull is
//  the 6DOF surface triangulation, the state that of the rigid-body core.
//
//  ship.dat (one keyword per line, # comments; all optional):
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
//  The ship frame is the body frame of the 6DOF object (x forward, z up at the start).
//  Output: REEF3D_<model>_6DOF/REEF3D_ship_<n>.dat

class ship : public sixdof_load
{
public:
    
    ship(lexer*, int);
    virtual ~ship();
    
    void add_load(lexer*, const sixdof_rigidbody&, const sixdof_geometry&, double*) override;
    void print(lexer*) override;
    
private:
    
    void read(lexer*);
    void ini(lexer*, const sixdof_rigidbody&, const sixdof_geometry&);
    
    const int id;
    bool initialized;
    
    // input
    double lpp, S, k, nu, B44, B44q, Cd, thrust, xthrust, zthrust;
    int friction, nstrip;
    bool lpp_in, S_in;
    
    // hull
    double xa, xf, zw;
    vector<double> xs, dx, T;
    
    // loads of the last evaluation, ship frame
    double ub, vb, wb, pb, qb, rb_;
    double Re, CF, XF, Ycf, Ncf, Kroll;
    
    ofstream out;
};

#endif
