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

/*--------------------------------------------------------------------
DEM core: rigid particles of arbitrary shape (level-set DEM) with
non-smooth contact dynamics (Moreau-Jean time stepping, velocity-level
impulses, Coulomb friction cone, projected Gauss-Seidel).

The core is independent of the flow solvers. Fluid coupling is done by
the adapter class dem_f, which sets the per-body external loads,
implicit drag coefficients and added mass before each step, and which
may supply extra contacts (e.g. against the REEF3D topo/solid level set).
--------------------------------------------------------------------*/

#ifndef DEM_CORE_H_
#define DEM_CORE_H_

#include<Eigen/Dense>
#include<vector>
#include<string>
#include<functional>
#include<unordered_map>
#include<cstdint>

using namespace std;

typedef Eigen::Vector3d dem_vec;
typedef Eigen::Matrix3d dem_mat;
typedef Eigen::Quaterniond dem_quat;

enum dem_shape_type
{
    DEM_SPHERE = 0,
    DEM_BOX = 1,
    DEM_CYLINDER = 2,
    DEM_ELLIPSOID = 3,
    DEM_MESH = 4
};

struct dem_material
{
    double rho = 2650.0;
    double friction = 0.5;
    double restitution = 0.3;
};

class dem_shape
{
public:
    dem_shape() = default;

    // builders, all shapes are centred at their centroid and aligned with their principal axes
    void build_sphere(double r);
    void build_box(double lx, double ly, double lz);
    void build_cylinder(double r, double length);             // axis along body z
    void build_ellipsoid(double a, double b, double c, int res);
    bool build_mesh(const string &file, double scale, int res); // ASCII or binary STL

    void make_nodes(double spacing);
    void make_quadrature(int n);

    // signed distance in body frame (negative inside) and its gradient
    double sdf(const dem_vec &y) const;
    dem_vec sdf_grad(const dem_vec &y) const;

    int type = DEM_SPHERE;
    dem_vec dim = dem_vec::Zero();      // sphere: (r,-,-) box: half extents, cylinder: (r,r,half length), ellipsoid: semi axes

    double volume = 0.0;
    double area = 0.0;
    double rbound = 0.0;                // bounding radius about centroid
    double deq = 0.0;                   // volume equivalent diameter
    double sphericity = 1.0;
    dem_vec inertia_unit = dem_vec::Zero();  // principal moments per unit density

    // surface nodes (body frame) used for level-set contacts
    vector<dem_vec> nodes;

    // volume quadrature (body frame) used for buoyancy
    vector<dem_vec> qp;
    vector<double> qw;

    // surface triangulation (body frame), used for output and meshes
    vector<dem_vec> vert;
    vector<int> tri;

private:
    void mesh_mass_properties(dem_vec &com, dem_mat &J, double &vol);
    void build_sdf_grid(int res);
    void finish_mesh_shape(bool rotate=true);
    double sdf_grid(const dem_vec &y) const;
    dem_vec sdf_grid_grad(const dem_vec &y) const;
    double mesh_distance(const dem_vec &y, double &winding) const;
    void icosphere(int levels);

    // SDF grid
    dem_vec g0 = dem_vec::Zero();
    double gdx = 1.0;
    int gn[3] = {0,0,0};
    vector<float> sdfval;
};

struct dem_body
{
    int shape = 0;
    int mat = 0;
    int id = 0;
    bool fixed = false;
    bool active = true;
    int mode = 0;               // 0 unresolved, 1 resolved (set by the adapter)

    double m = 0.0;             // mass
    dem_vec Ib = dem_vec::Zero(); // principal inertia (body frame)

    dem_vec x = dem_vec::Zero();
    dem_quat q = dem_quat::Identity();
    dem_vec v = dem_vec::Zero();
    dem_vec w = dem_vec::Zero(); // angular velocity, world frame
    dem_vec vold = dem_vec::Zero();
    dem_vec wold = dem_vec::Zero();
    dem_vec vp = dem_vec::Zero();   // split-impulse pseudo velocities (position correction only)
    dem_vec wp = dem_vec::Zero();

    // loads set by the adapter before each step (not including gravity)
    dem_vec F = dem_vec::Zero();
    dem_vec T = dem_vec::Zero();

    // implicit drag: F_d = K (uf - v), T_d = Kr (-w), added mass madd with fluid acceleration af
    double K = 0.0;
    double Kr = 0.0;
    double madd = 0.0;
    dem_vec uf = dem_vec::Zero();
    dem_vec af = dem_vec::Zero();

    // work
    dem_vec invm = dem_vec::Zero();  // diagonal inverse mass in world frame (2D masks)
    dem_mat invI = dem_mat::Zero();  // inverse inertia in world frame
    dem_mat R = dem_mat::Identity();
};

struct dem_contact
{
    int a = -1;                 // body receiving +P
    int b = -1;                 // second body, <0 for a wall
    uint64_t key = 0;
    dem_vec x = dem_vec::Zero();
    dem_vec n = dem_vec::UnitZ();  // unit normal pointing from b to a
    double gap = 0.0;
    double mu = 0.5;
    double e = 0.3;
    dem_vec vwall = dem_vec::Zero();

    // solver work
    dem_vec t1,t2,ra,rb;
    dem_mat Wl;
    dem_vec P = dem_vec::Zero();    // local impulse (n,t1,t2)
    double target = 0.0;
    double stab = 0.0;
    double Pp = 0.0;            // split-impulse normal pseudo impulse
    double un0 = 0.0;
    bool speculative = false;
};

struct dem_plane
{
    dem_vec n = dem_vec::UnitZ();   // unit normal pointing into the domain
    double d = 0.0;                 // plane: n.x = d
};

class dem_core
{
public:
    dem_core() = default;

    // wall callback: adds contacts for the given margin; must give identical results on all ranks
    typedef std::function<void(dem_core&, double, vector<dem_contact>&)> wallfunc;

    void initialize();
    void step(double dt, const wallfunc &walls);

    double kinetic_energy() const;
    double max_velocity() const;
    double min_rbound() const;

    uint64_t make_key(int a, int b, int feature) const;

    vector<dem_shape> shapes;
    vector<dem_material> mats;
    vector<dem_body> bodies;
    vector<dem_plane> planes;
    dem_material wallmat;

    dem_vec gravity = dem_vec(0.0,0.0,-9.81);
    bool plane2D = false;       // restrict motion to the x-z plane (REEF3D 2D: j_dir==0)

    // solver parameters
    int maxiter = 1000;
    double tol = 1.0e-5;
    double theta = 0.5;         // Moreau-Jean midpoint for contact activation
    double beta = 0.2;          // penetration stabilisation
    double slop = 1.0e-5;       // allowed penetration [m]
    double vstab_max = 0.5;     // max stabilisation velocity [m/s]
    double vrest = 1.0e-3;      // restitution threshold velocity
    double warmstart = 1.0;
    int manifold = 6;           // max contacts kept per body pair and normal cluster (0: keep all)
    int pseudoiter = 50;        // split-impulse iterations
    double sor = 1.0;           // relaxation factor of the Gauss-Seidel updates

    // statistics of the last step
    int ncontacts = 0;
    int iterations = 0;
    double residual = 0.0;
    double maxpen = 0.0;

    vector<dem_contact> contacts;

    // helpers for wall/other contact generators
    dem_vec node_world(const dem_body &B, const dem_vec &y) const {return B.x + B.R*y;}

private:
    void update_mass(double dt);
    void detect(double margin);
    void narrow(int a, int b, double margin);
    void planes_contacts(double margin);
    void reduce();
    void prepare(double dt);
    void solve(double dt);
    void integrate(double dt);
    dem_vec gyro_free(const dem_body &B, double dt) const;

    void apply(int c, const dem_vec &Pw);
    void apply_pseudo(int c, double Pn);
    dem_vec relvel(const dem_contact &C) const;

    unordered_map<uint64_t,dem_vec> cacheP;     // world impulse of last step
    unordered_map<uint64_t,double> cacheU;      // approach velocity of speculative contacts
};

#endif
