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
// REEF3D::CORAL: soil-structure interaction of piles with nonlinear springs (Winkler foundation).
//
// soil_curve: backbone of a spring as an odd piecewise-linear table (y >= 0; force for -y is -force(y);
// flat beyond the last point), or compression only (Q-z). API RP 2GEO / RP 2A curves:
//   p-y sand (Reese, Cox & Koop 1974 as in API): p = A pu tanh(k X y/(A pu)),
//       A = max(0.9, 3 - 0.8 X/D) (static) or 0.9 (cyclic), pu = min(p_wedge, p_flow) per length
//   p-y soft clay (Matlock 1970): pu = min(3 su + gamma' X + J su X/D, 9 su) D, yc = 2.5 eps50 D,
//       p/pu = 0.5 (y/yc)^(1/3) up to 8 yc (static); cyclic: 0.72 pu beyond 3 yc below XR, degrading to
//       0.72 pu X/XR at 15 yc above XR = 6 D/(gamma' D/su + J); the infinite initial slope is
//       replaced by the secant to 0.01 yc
//   t-z clay: t_max = alpha su, alpha = 0.5 psi^-0.5 (psi <= 1) or 0.5 psi^-0.25 (psi = su/p0'),
//       t/t_max at z/D = 0.0016 0.0031 0.0057 0.0080 0.0100 0.0200: 0.3 0.5 0.75 0.9 1.0 residual
//   t-z sand: t_max = K p0' tan(delta), full at z = 2.54 mm
//   Q-z: q = 9 su (clay) or Nq p0' (sand), Q/Qp at z/D = 0.002 0.013 0.042 0.073 0.100:
//       0.25 0.5 0.75 0.9 1.0 (compression only)
// soil_profile: layers below the mudline (sand: phi, gamma', k, delta, K, Nq; clay: su linear in the
// layer, gamma', eps50, J; linear: Winkler moduli kh, kv and tip stiffness kq), effective stress from the
// layers above.
// soil_spring: the springs of one pile node (a frame node body): lateral p-y (radial, isotropic about
// the pile axis), axial t-z and tip Q-z along the axis, each times the node's tributary length, and the
// torsion of the shaft friction (torque t(r phi) r, r = D/2, phi the twist about the axis); displacements
// from the node's reference position. Nonlinear elastic (no hysteresis); analytic tangent.

#ifndef CORAL_SOIL_H_
#define CORAL_SOIL_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_types.h"
#include<string>
#include<vector>

namespace coral
{
    class model;

    struct soil_curve
    {
        std::vector<double> y, p;      // y[0] = 0, p[0] = 0, increasing y
        bool compression_only = false; // force only for y > 0 (Q-z: y = settlement)

        double force(double y) const;
        double slope(double y) const;
        bool empty() const {return y.size() < 2;}
        soil_curve scaled(double s) const;          // force times s
        static soil_curve linear(double k, double ymax=1.0e3);
    };

    struct soil_layer
    {
        enum class kind {sand, clay, linear};
        kind type = kind::sand;
        double top = 0.0, bottom = 0.0;             // depths below the mudline [m]
        double gamma = 10.0e3;                      // submerged unit weight [N/m^3]
        // sand
        double phi = 30.0;                          // friction angle [deg]
        double k = 0.0;                             // initial modulus of subgrade reaction [N/m^3]
        double delta = -1.0;                        // pile-soil friction angle [deg] (< 0: phi - 5)
        double K = 0.8;                             // lateral earth pressure coefficient (t-z)
        double Nq = 40.0;                           // bearing capacity factor
        // clay
        double su_top = 0.0, su_bottom = -1.0;      // undrained shear strength [Pa] (bottom < 0: constant)
        double eps50 = 0.01;
        double J = 0.5;
        double residual = 0.9;                      // t-z residual t/t_max
        // linear (Winkler): lateral and axial modulus per length [N/m^2], tip stiffness [N/m]
        double kh = 0.0, kv = 0.0, kq = 0.0;
        bool cyclic = false;

        double su(double depth) const;
    };

    class soil_profile
    {
    public:
        std::string name;
        double mudline = 0.0;                       // z of the mudline
        std::vector<soil_layer> layers;

        const soil_layer *layer_at(double depth) const;
        double effective_stress(double depth) const;     // p0' [Pa]
        soil_curve py(double depth, double D) const;     // lateral force per length against y
        soil_curve tz(double depth, double D) const;     // axial force per length against z
        soil_curve qz(double depth, double D) const;     // tip force against settlement (closed tip)
    };

    class soil_spring : public load
    {
    public:
        // node body b at its current position, pile axis a; lateral and axial curves are already per
        // node (times the tributary length), tip curve only at the pile tip (else empty)
        soil_spring(const model&, int b, const Vec3 &axis, const soil_curve &lateral, const soil_curve &axial,
                    const soil_curve &tip, const std::string &name, double radius=0.0);

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
        void eval(const model&, double t, VecX &f) const override;
        bool has_tangent() const override {return true;}
        void tangent(const model&, double t, double cK, double cC, double cM, std::vector<triplet>&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // soil force on the node (global)
        Vec3 force(const model&) const;

    private:
        std::string name_;
        int b_;
        Vec3 a_, x0_;
        Mat3 R0_;
        double r_;
        soil_curve lat_, ax_, tip_;
    };

    // springs along a pile made of frame nodes (bodies at the given stations, top to tip or tip to top):
    // nodes below the mudline get p-y (and t-z, Q-z at the deepest node) with tributary lengths;
    // linear: the initial slopes only. Returns the number of springs added.
    int add_pile_springs(model&, const soil_profile&, const std::vector<int> &nodes, double D, bool py, bool tz,
                         bool qz, bool linear, const std::string &name);
}

#endif
