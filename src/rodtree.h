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

#ifndef RODTREE_H_
#define RODTREE_H_

// Branching flexible rods (soft corals, gorgonians, flexible vegetation).
//
// Discrete Cosserat rod on a tree: every rod element is a rigid cylinder
// (centre c, orientation q, body z-axis = rod axis) and neighbouring
// elements are connected by elastic 6-DOF joints at the shared node:
//   translation: penalty spring  k = E A / l_j        (extension/shear)
//   rotation   : K = EI/l_j (I - t t^T) + GJ/l_j t t^T acting on the
//                rotation vector of the relative rotation w.r.t. rest
// Kelvin-Voigt damping (time constant beta) on both.  Branch points are
// ordinary joints between the parent element and each child element, so
// bending, torsion and junction stiffness follow from the same model.
// Root nodes are clamped to the bed (joint to a fixed world frame).
//
// Hydrodynamics (unresolved, Morison-type, per element with relative
// velocity):  normal/tangential drag, Froude-Krylov + added-mass inertia
// Cm*rho*V*a_f,n, added mass Ca*rho*V in the normal directions of the
// element mass matrix, buoyancy, gravity.  All fluid terms scale with the
// submerged fraction chi supplied by the coupling.
//
// Time integration: linearly implicit Euler on the full tree (sparse LU,
// joint Jacobians by central finite differences, drag linearised
// analytically) or symplectic Euler with automatic sub-cycling.
//
// This file has no REEF3D dependencies (Eigen + STL only) so the solver
// can be tested standalone (tests/rodtree).

#include<vector>
#include<string>
#include<iosfwd>
#include<Eigen/Dense>
#include<Eigen/Sparse>

class rodtree
{
public:
    typedef Eigen::Vector3d Vec3;
    typedef Eigen::Matrix3d Mat3;
    typedef Eigen::Quaterniond Quat;

    struct material
    {
        double E=1.0e7, nu=0.4, rho=1100.0;    // Pa, -, kg/m3
        double Cdn=1.2, Cdt=0.01;                // normal / tangential drag
        double Ca=1.0, Cm=2.0;                   // added mass, inertia coeff.
        double beta=0.0;                         // Kelvin-Voigt time constant [s]
        // polyps: porous layer of extended polyps around the branch
        double hp=0.0, phip=0.0;                 // polyp height [m], frontal solidity of the layer [-]
        double Cdp=1.0, ktp=1.0;                 // polyp drag coeff., axial/normal ratio of polyp drag
        double Ur=-1.0, dUr=0.0, taup=0.0;       // retraction: onset speed, ramp width [m/s], time constant [s]
    };

    struct element
    {
        int colony, n0, n1;                      // node ids (input numbering, after refinement)
        double l, r;                             // length, radius
        double m;                                // structural mass
        Vec3 c, v, w;                            // centre, velocity, angular velocity (world)
        Quat q;                                  // orientation, body z = axis
        Vec3 c0; Quat q0;                        // initial configuration
        Vec3 a;                                  // last structural acceleration (for reaction)
        // fluid state supplied by the coupling
        Vec3 uf, af; double chi;
        // hydrodynamic force of the last load evaluation (on the structure)
        Vec3 Fh, Fdrag, Finert;
        // polyp extension 0 (retracted) .. 1 (fully extended)
        double ext;
    };

    struct joint
    {
        int a, b;                                // parent element (-1 = world), child element
        Vec3 pa, pb;                             // attachment in body frames (pa world point if a<0)
        Quat qrel0;                              // rest relative rotation conj(qa)*qb
        double k, kb, kt;                        // translational, bending, torsional stiffness
        double beta;
        Vec3 F, T;                               // last force/torque on child (reaction at root)
    };

    struct colony_info
    {
        std::string name;
        material mat;
        std::vector<int> elements, roots;        // element ids, root joints
        int tip_node=-1;                         // element whose end is farthest from root
    };

    rodtree();

    // input
    void read(std::istream&);                    // throws std::runtime_error on bad input
    void finalize_setup();                       // build elements/joints from parsed data

    // simulation
    void set_gravity(const Vec3& g_){g=g_;}
    void set_fluid_density(double r){rhof=r;}
    void set_integrator(int m){integrator=m;}   // 0 implicit, 1 explicit
    void set_substeps(int n){substeps=n>0?n:1;}
    int get_integrator() const {return integrator;}
    int get_substeps() const {return substeps;}
    int get_reaction_mode() const {return reaction_mode;}  // 0 full, 1 drag only, 2 one-way
    void compute_hydro();                        // Fh, Fdrag, Finert with current state
    void advance(double dt);                     // one fluid time step

    // reaction force on the fluid of element e (Newton); by reaction_mode:
    // full = -(drag + inertia), drag = -drag, none = one-way coupling
    Vec3 fluid_reaction(int e) const;
    // d|F|/d|u_rel| of the element drag (for point-implicit fluid drag)
    double drag_slope(int e) const;

    // queries
    int nelem() const {return (int)el.size();}
    int njoint() const {return (int)jt.size();}
    int ncolony() const {return (int)col.size();}
    const element& elem(int e) const {return el[e];}
    element& elem(int e) {return el[e];}
    const joint& jnt(int j) const {return jt[j];}
    const colony_info& colony(int c) const {return col[c];}
    Vec3 axis(int e) const {return el[e].q*Vec3(0,0,1);}
    Vec3 end0(int e) const {return el[e].c - 0.5*el[e].l*axis(e);}
    Vec3 end1(int e) const {return el[e].c + 0.5*el[e].l*axis(e);}
    Vec3 point(int e, int iq, int nq) const {return el[e].c + (-0.5 + (iq+0.5)/double(nq))*el[e].l*axis(e);}
    double volume(int e) const;
    double time() const {return t;}
    int substeps_used() const {return nsub_last;}

    // diagnostics
    double kinetic_energy() const;
    double elastic_energy() const;
    Vec3 tip_displacement(int c) const;
    Vec3 base_force(int c) const;                // sum of root joint forces on the structure
    Vec3 hydro_force(int c) const;
    double polyp_extension(int c) const;         // mean over the colony's elements

    // output
    void write_vtp(const std::string& filename) const;

private:
    // parsed input
    struct pnode {int id; Vec3 x; double r;};
    struct pedge {int p, c;};
    struct pcolony
    {
        std::string name; material mat; int refine=1;
        std::vector<pnode> nodes; std::vector<pedge> edges; std::vector<int> clamps;
        std::vector<Vec3> instances;
    };
    std::vector<pcolony> input;

    std::vector<element> el;
    std::vector<joint> jt;
    std::vector<colony_info> col;

    Vec3 g;
    double rhof, t;
    int integrator, substeps, nsub_last, reaction_mode;
    double hsub_prev;

    // implicit solver data
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > lu;
    bool pattern_ready;
    std::vector<Eigen::Triplet<double> > trip;

    // helpers
    static Vec3 logq(const Quat&);
    static Quat expq(const Vec3&);
    Mat3 inertia_world(int e, bool with_fluid) const;
    Mat3 mass_trans(int e) const;
    void joint_force(const joint&, const Vec3& ca, const Quat& qa, const Vec3& va, const Vec3& wa,
                     const Vec3& cb, const Quat& qb, const Vec3& vb, const Vec3& wb,
                     Vec3& Fa, Vec3& Ta, Vec3& Fb, Vec3& Tb) const;
    void assemble_forces(std::vector<Vec3>& F, std::vector<Vec3>& T, bool store);
    void external_loads(int e, Vec3& F, Vec3& T, Mat3& Cdrag) const;
    void drag_coeffs(int e, double& cn, double& ct) const;
    double polyp_target(int e) const;
    void update_polyps(double dt);
    void step_implicit(double h);
    void step_explicit(double h);
    double explicit_dt() const;
};

#endif
