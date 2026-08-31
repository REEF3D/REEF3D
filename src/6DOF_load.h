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

#ifndef SIXDOF_LOAD_H_
#define SIXDOF_LOAD_H_

#include<vector>
#include<Eigen/Dense>

class lexer;
class sixdof_rigidbody;
class sixdof_geometry;

//  Fluid access of a load model, provided by the solver coupling (nullptr if the coupling has none):
//   velocity: fluid velocity at n points xyz[3n] (inertial frame); all ranks call it with the
//             same points and get the same result (points outside the domain: 0)
class sixdof_fluid
{
public:
    virtual ~sixdof_fluid() = default;
    
    virtual void velocity(int n, const double *xyz, double *uvw)=0;
};

//  Actuator disk (body-force propeller) as a momentum source of the fluid, inertial frame.
//  The propeller pushes the fluid with the force -T along the axis and swirls it with the
//  torque Q in the direction of the blade motion; the body gets +T along the axis (and the
//  shaft torque reaction, applied by the load model itself). The coupling distributes force
//  and torque over the cells inside the disk with the Hough-Ordway radial distribution,
//  normalised on its own grid so that the discrete totals are exactly T and Q.
struct sixdof_actuator_disk
{
    Eigen::Vector3d centre;     // disk centre
    Eigen::Vector3d axis;       // unit vector, direction of the thrust on the body
    double R, Rh;               // tip and hub radius
    double thickness;           // axial extent of the source region
    double T, Q;                // thrust [N] and torque [Nm] of the propeller
    int sense;                  // +1: blades turn right-handed about the axis, -1: left-handed
    
    // Hough-Ordway weights at point x: axial weight wa >= 0 (thrust density ~ wa along the axis),
    // tangential weight wt >= 0 (swirl density ~ wt along et) and et, the unit vector of the blade
    // motion, with the distance r from the axis; false outside the disk
    bool weights(const Eigen::Vector3d &x, double &wa, double &wt, Eigen::Vector3d &et, double &r) const;
    
    // discrete distribution on the grid of the coupling, per force component c (the source points
    // of component c, e.g. the staggered velocity points of CFD or the cell centres of NHFLOW):
    //  accumulate: sums over the points x with volume V and fluid fraction fl (0..1) of the cell;
    //              the coupling sums them over the MPI ranks
    //  swirl_factor: scale of the tangential force from the sums of all three components
    //  force:      force density of component c [N/m^3] at point x
    // The thrust part gives exactly -T along the axis on each grid. The swirl part has no net
    // force on each grid (the tangential unit vectors of the few cells of a disk on a coarse grid
    // do not cancel by themselves), and with the thrust part (which has a small torque when the
    // components sit on different grids) exactly the torque Q about the axis.
    struct sums
    {
        double SA0=0.0, SA=0.0, SW=0.0, SE=0.0, G1=0.0, G0=0.0, GA=0.0;
    };
    
    void accumulate(const Eigen::Vector3d &x, double V, double fl, int c, sums&) const;
    double swirl_factor(const sums*) const;
    double force(const Eigen::Vector3d &x, double fl, int c, const sums&, double kappa) const;
};

//  External load model acting on a 6DOF body (solver independent), e.g. the ship module:
//  evaluated with the body state of every stage, before the rigid-body right-hand side.
//
//  add_load: add the load to F[0..5] = X, Y, Z, K, M, N (inertial frame, moments about the
//            centre of gravity); fluid: velocity sampling of the coupling (may be nullptr)
//  actuator_disks: momentum sources of the fluid (body-force propellers) of the last add_load
//  added_mass: adds the added mass of the model to A (body frame, about the centre of gravity,
//            order surge, sway, heave, roll, pitch, yaw); true if the model has one. The coupling
//            solves the equations of motion with M + A (implicit, as for the FNPF added mass).
//  fluid_mask: weights w[6] of the hydrodynamic loads of the solver coupling (inertial X, Y, Z,
//            K, M, N); the model sets entries to 0 when it supplies these loads itself
//  print:    output once per time step (rank 0 decides inside)

class sixdof_load
{
public:
    virtual ~sixdof_load() = default;
    
    virtual void add_load(lexer*, const sixdof_rigidbody&, const sixdof_geometry&, sixdof_fluid*, double*)=0;
    virtual void actuator_disks(std::vector<sixdof_actuator_disk>&) const {}
    virtual bool added_mass(const sixdof_rigidbody&, Eigen::Matrix<double,6,6>&) const {return false;}
    virtual void fluid_mask(double*) const {}
    virtual void print(lexer*) {}
};

#endif
