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

#ifndef ICE_CONTACT_NSCD_H_
#define ICE_CONTACT_NSCD_H_

#include"ice_contact.h"
#include<map>
#include<tuple>

using namespace std;

// Non-smooth contact dynamics for convex polygons in the horizontal plane (Moreau-Jean type):
// velocity-level Signorini condition with Newton restitution and Coulomb friction, solved per
// time step by projected Gauss-Seidel on the contact impulses (sequential impulses), warm started
// from the previous step. Contacts are speculative: a gap g >= 0 closing faster than g/dt is
// stopped, so large time steps do not tunnel. Penetration is removed by a Baumgarte bias.
//
// Narrow phase: separating axis test, reference face / incident edge clipping, up to two
// contact points per polygon pair. Walls: the domain box as fixed half-planes (optional).

class ice_contact_nscd final : public ice_contact
{
public:
    ice_contact_nscd(double mu, double e, int iter, int walls, double xmin, double xmax, double ymin, double ymax);
    virtual ~ice_contact_nscd() = default;

    void solve(vector<ice_body2D>&, double dt, int is2D) override final;
    const vector<ice_contact_record>& records() const override final {return rec;}

private:
    struct contact
    {
        int a,b;              // bodies, b = -1-wall for a wall
        double nx,ny;         // normal, from a to b
        double px,py;         // contact point
        double sep;           // signed gap (<0 penetration)
        long key;             // feature key for warm starting
        double rax,ray,rbx,rby;
        double mn,mt;         // effective masses
        double bias;          // target normal velocity
        double Pn = 0.0, Pt = 0.0;
    };

    void broadphase(vector<ice_body2D>&, double dt);
    void collide(vector<ice_body2D>&, int, int, double);
    void walls(vector<ice_body2D>&, int, double);
    double max_separation(const ice_body2D&, const ice_body2D&, int&) const;

    double mu,e;
    int iter;
    int wallflag;
    double xmin,xmax,ymin,ymax;
    double beta,slop,vrest;

    vector<contact> con;
    vector<ice_contact_record> rec;
    map<tuple<int,int,long>,pair<double,double>> cache;
};

#endif
