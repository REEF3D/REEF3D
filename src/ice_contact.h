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

#ifndef ICE_CONTACT_H_
#define ICE_CONTACT_H_

#include<vector>

using namespace std;

// Planar rigid body for the contact solver: convex polygon in world coordinates,
// velocities at the centre of mass. im, iI = 0 marks a fixed body.
// Kept free of the FNPF ice data so a standalone REEF3D DEM can take this slot later.

struct ice_body2D
{
    double x = 0.0, y = 0.0;            // centre of mass
    double vx = 0.0, vy = 0.0, w = 0.0; // velocity, yaw rate
    double im = 0.0, iI = 0.0;          // inverse mass, inverse yaw inertia
    double rb = 0.0;                    // bounding radius about (x,y)
    vector<double> px,py;               // convex polygon, counter-clockwise, world
    double fx = 0.0, fy = 0.0;          // out: contact force on the body (impulse/dt)
    int ncontact = 0;                   // out
    int id = 0;                         // persistent id, keys the warm start
};

// Active contact after a solve, per body pair (b = -1-wall for walls): contact point, normal from a to b,
// total normal force (normal impulse/dt, summed over the pair's contact points).
struct ice_contact_record
{
    int a = 0, b = 0;
    double px = 0.0, py = 0.0;
    double nx = 0.0, ny = 0.0;
    double Fn = 0.0;
};

class ice_contact
{
public:
    virtual ~ice_contact() = default;
    virtual void solve(vector<ice_body2D>&, double dt, int is2D)=0;
    virtual const vector<ice_contact_record>& records() const=0;
};

#endif
