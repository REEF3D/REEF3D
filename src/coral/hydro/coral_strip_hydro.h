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
// REEF3D::CORAL: strip-theory hydrodynamics of a rigid floating body built from cylinders (columns, pontoons, braces of
// a semi-submersible or a spar), with nonlinear hydrostatics: the engineering model of a floater without a
// hydrodynamic database, as OpenFAST's strip theory.
// Per strip of a member (length L, diameter D, A = pi D^2/4), with the wet fraction f of the strip below the
// instantaneous wave elevation (continuous in time):
//   buoyancy   rho g A f L, upwards, at the centre of the wet part (the restoring of heave, roll and pitch follows from
//              the changing wet length of the surface-piercing members and the moving centres)
//   Morison    f L (rho A Cm a_w,n - rho A Ca a_b,n + 1/2 rho Cd D |u_r,n| u_r,n), normal to the axis, u_r = u_w + U_c - u_b
// At a closed member end below the surface (e.g. the bottom of a column; flag end1/end2), along the axis:
//   pressure   -rho p_dyn/rho A n_out (the dynamic pressure of the waves on the end face: heave excitation)
//   added mass Ca_end rho D^3/6 (a_w,ax - a_b,ax) (Ca_end 1: half the added mass of a disk, rho D^3/3)
//   drag       1/2 rho Cd_end A |u_r,ax| u_r,ax
// Wave kinematics from the wave source (CORAL's linear waves or a sampled fluid) at the strip centres and the end
// points where they were at the start, once per time (as OpenFAST's HydroDyn: the undisplaced positions; the floater's
// motion is small against the wave length), or at their current positions with exact_kinematics (each evaluation, as
// the time integrator perturbs and iterates the state: costly with many wave components); above the still water level
// those of z = 0.

#ifndef CORAL_STRIP_HYDRO_H_
#define CORAL_STRIP_HYDRO_H_

#include"coral/core/coral_load.h"
#include"coral/hydro/coral_waves.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    struct strip_member
    {
        Vec3 p1 = Vec3::Zero(), p2 = Vec3::Zero();    // global, at creation
        double D = 1.0, Cd = 1.0, Ca = 1.0, Cm = 2.0;
        int n = 20;
        bool end1 = false, end2 = false;              // closed end faces with pressure, added mass and drag
        double Ca_end = 1.0, Cd_end = 1.0;
        bool buoyancy = true;
    };

    class strip_hydro : public load
    {
    public:
        strip_hydro(const model&, int body, const std::vector<strip_member>&, std::shared_ptr<const wave_kinematics>,
                    double rho, double g, const Vec3 &current, const std::string &name);

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
        void eval(const model&, double t, VecX &f) const override;
        void kinematics_points(const model&, const void *source, std::vector<Vec3>&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        bool exact_kinematics = false;

        // total force (global) and moment about the body's CoG (global), submerged volume, of the current state
        void total(const model&, double t, Vec3 &F, Vec3 &M, double &V) const;

    private:
        struct strip {Vec3 s; Vec3 e; double L, D, Cd, Ca, Cm; bool buoy;};
        struct end {Vec3 s; Vec3 n; double D, Ca, Cd;};     // body frame: point, outward normal
        // the forces of the current state at their points (callback per force)
        template<class Add> void forces(const model&, double t, Add &&add, double &V) const;

        std::string name_;
        int b_;
        std::vector<strip> strips_;
        std::vector<end> ends_;
        std::shared_ptr<const wave_kinematics> waves_;
        double rho_, g_;
        Vec3 current_;
        // the kinematics at the undisplaced points (strips, then ends), cached for one time
        struct kin {double eta; Vec3 u, du; double p;};
        std::vector<Vec3> p0_;
        mutable std::vector<kin> cache_;
        mutable double t_cache_ = -1.0e300;
        void kinematics(int idx, const Vec3 &p, double t, kin &k) const;
    };
}

#endif
