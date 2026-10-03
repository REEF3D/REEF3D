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

#ifndef FEM_SOLID_H_
#define FEM_SOLID_H_

// Explicit Total-Lagrangian finite element solver for solids
// (elastic deformation, stresses, material failure and collapse).
//
//  - mesh: 8-node hexahedra on a uniform voxel lattice, filled from boxes
//    and closed STL surfaces (shared lattice nodes glue touching bodies)
//  - element: hex8, full 2x2x2 integration or one-point integration with
//    Flanagan-Belytschko stiffness hourglass control
//  - materials: elastic (St. Venant-Kirchhoff), J2 plasticity with linear
//    hardening and failure strain, concrete (Rankine tension / compression
//    crushing isotropic damage, crack-band regularised by G_f and G_c)
//  - failure: element erosion, the mass stays on the nodes; nodes without
//    any intact element become discrete debris particles
//  - contact: ground plane and node-node penalty contact between parts
//    that do not share an intact element (fragments, debris, separate bodies)
//  - time integration: central differences with lumped mass, the solver
//    subcycles inside a given (fluid) time step
//
// The class has no REEF3D or MPI dependencies (Eigen only). The coupling
// to the flow solver is in fem_coupling.

#include<vector>
#include<string>
#include<iosfwd>
#include<algorithm>
#include<Eigen/Dense>

class fem_solid
{
public:

    typedef Eigen::Vector3d Vec3;
    typedef Eigen::Matrix3d Mat3;

    enum {MAT_ELASTIC=0, MAT_J2=1, MAT_CONCRETE=2};

    struct material
    {
        int id = 0;
        int type = MAT_ELASTIC;
        double rho = 1000.0, E = 1.0e6, nu = 0.3;
        double lambda = 0.0, mu = 0.0, cp = 0.0;
        // J2
        double sigy = 0.0, H = 0.0, epsfail = 0.0;
        // concrete
        double ft = 0.0, Gf = 0.0, fc = 0.0, Gc = 0.0;
        // erosion threshold for damage
        double derode = 0.99;
    };

    // history variables of one integration point
    struct gpstate
    {
        double kt = 0.0, kc = 0.0;      // damage history (equivalent strains)
        double d = 0.0;                 // total damage
        double ep = 0.0;                // equivalent plastic strain
        double Ep[6] = {0,0,0,0,0,0};   // plastic Green-Lagrange strain (xx yy zz xy yz zx)
    };

    struct element
    {
        int n[8];
        int mat = 0;                    // index into mats
        int ix, iy, iz;                 // voxel index
        bool alive = true;
        double svm = 0.0;               // von Mises (Cauchy) stress, averaged over the GPs
        double J = 1.0;
    };

    // surface face of an intact element, nodes ordered counter-clockwise seen
    // from outside (normal = (x2-x0) x (x3-x1) points outward)
    struct face
    {
        int n[4];
        int elem;
    };

    struct monitor
    {
        std::string name;
        int node;
    };

    // coupling options (read from the same input file, used by fem_coupling)
    struct coupling_options
    {
        double pressure_offset = 1.5;   // pressure probe distance from the surface [fluid cells]
        int    points_per_face = 0;     // Lagrangian points per face edge, 0: automatic
        double debris_cd = 1.0;         // drag coefficient of debris particles
        int    debris_reaction = 1;     // 1: debris drag acts back on the fluid
        int    forcing = 1;             // 1: direct forcing of the surface velocity on the fluid
        double print_dt = 0.0;          // VTU print interval, overrides Z 31 if > 0
        int    loads = 0;               // 0: implicit direct-forcing reaction, 1: pressure integration (explicit)
    };

    fem_solid();
    ~fem_solid();

    // ------------------------------------------------------------------
    // setup
    // ------------------------------------------------------------------
    void read(std::istream&);                       // fem.dat (throws std::runtime_error)
    void build();                                   // voxelise, nodes, mass, surface
    void set_gravity(const Vec3& g) {grav = g;}
    void set_plane_strain(bool b) {plane_strain = b;}

    // programmatic setup (tests): same as the input keywords
    void set_lattice(double ox,double oy,double oz,double hx,double hy,double hz);
    int  add_material(const material&);
    void add_box(double x0,double x1,double y0,double y1,double z0,double z1,int matid);
    void remove_box(double x0,double x1,double y0,double y1,double z0,double z1);
    void add_fix(double x0,double x1,double y0,double y1,double z0,double z1,bool fx,bool fy,bool fz);
    void set_element_type(int full) {full_int = full;}
    void set_hourglass(double c) {hg_coef = c;}
    void set_cfl(double c) {cfl = c;}
    void set_damping(double a) {alpha_damp = a;}
    void set_ground(double z,double kfac,double mu) {ground_on = true; zground = z; kground = kfac; mu_ground = mu;}
    void set_contact(bool on,double kfac,double mu) {contact_on = on; kcontact = kfac; mu_contact = mu;}

    // ------------------------------------------------------------------
    // time stepping
    // ------------------------------------------------------------------
    // advance by dt with automatic subcycling, external (fluid) loads are
    // kept constant over the step
    void advance(double dt);
    int  last_substeps() const {return nsub_last;}
    double critical_dt() const {return dtcrit;}
    double time() const {return t;}

    // ------------------------------------------------------------------
    // state access (coupling)
    // ------------------------------------------------------------------
    int nnode() const {return (int)X.size();}
    int nelem() const {return (int)elems.size();}
    const Vec3& ref_pos(int i) const {return X[i];}
    const Vec3& pos(int i) const {return x[i];}
    const Vec3& vel(int i) const {return v[i];}
    double mass(int i) const {return m[i];}
    double node_volume(int i) const {return vnode[i];}
    double node_density(int i) const {return vnode[i]>0.0 ? m[i]/vnode[i] : 0.0;}
    void set_vel(int i,const Vec3& vv) {v[i] = vv;}
    void set_pos(int i,const Vec3& xx) {x[i] = xx;}
    const element& elem(int e) const {return elems[e];}
    const gpstate& gp(int e,int q) const {return gps[e*ngp+q];}
    int  ngauss() const {return ngp;}
    const std::vector<face>& surface() const {return faces;}
    const std::vector<int>& debris() const {return orphan;}
    bool is_debris(int i) const {return nalive[i]==0;}
    int  surface_version() const {return surf_version;}
    double hmin() const {return std::min(hx,std::min(hy,hz));}
    double lattice_h(int d) const {return d==0 ? hx : d==1 ? hy : hz;}
    const coupling_options& coupling() const {return copt;}
    const std::vector<monitor>& monitors() const {return mons;}

    void clear_loads();
    void add_load(int i,const Vec3& f) {fext[i] += f;}

    // fluid coupling, set by the coupling once per fluid step:
    //  - the fluid parcel in the forcing volume of the surface points (mass M,
    //    momentum M u_f) is merged with the node at the start of the step and
    //    moves with it during the substeps: (m - m_f + M) a = F + (m - m_f) g
    //  - m_f: fluid mass in the node volume (buoyancy and inertia of the fluid
    //    enclosed by the immersed boundary)
    void clear_coupling();
    void add_coupling(int i,double M,const Vec3& Mu) {Mcpl[i] += M; Mu_cpl[i] += Mu;}
    void set_fluid_mass(int i,double mf) {mfl[i] = mf;}
    double extent() const;                      // largest dimension of the initial geometry

    // ------------------------------------------------------------------
    // diagnostics
    // ------------------------------------------------------------------
    int n_alive() const;
    int n_eroded() const {return nelem()-n_alive();}
    double kinetic_energy() const;
    double strain_energy() const;               // elastic strain energy (undamaged part)
    double dissipated_energy() const {return wdiss;}
    Vec3 support_force() const;                 // force of the structure on its supports
    Vec3 support_force(double x0,double x1,double y0,double y1,double z0,double z1) const;   // supports inside a box
    Vec3 total_load() const;                    // sum of the external (fluid) loads
    double max_vonmises() const;
    double max_displacement() const;

    // ------------------------------------------------------------------
    // output
    // ------------------------------------------------------------------
    void write_vtu(const std::string& filename) const;
    void info(std::ostream&) const;

private:

    // mesh
    void voxelise();
    void voxelise_stl(const std::string& file,int matid,bool remove);
    void make_nodes_elements();
    void build_surface();
    void compute_mass();
    void count_bodies();
    int  voxel(int ix,int iy,int iz) const;

    // element / material
    void shape_derivatives();
    void internal_forces(double dts);
    void stress(const material&,gpstate&,const Mat3& F,const Mat3& Fdot,double h,double w,Mat3& P,double& svm,bool& failed);
    double damage_exp(double kappa,double e0,double ef) const;

    // contact
    void contact_ground();
    void contact_nodes();
    bool share_alive_element(int a,int b) const;

    // lattice
    double ox=0.0, oy=0.0, oz=0.0, hx=0.0, hy=0.0, hz=0.0;
    int nx=0, ny=0, nz=0;
    std::vector<int> vox;                       // voxel -> material index (-1 empty)
    std::vector<int> vox_elem;                  // voxel -> element index (-1 empty)

    struct shape_box {double x0,x1,y0,y1,z0,z1; int mat; bool remove;};
    struct shape_stl {std::string file; int mat; bool remove;};
    struct fix_box {double x0,x1,y0,y1,z0,z1; bool f[3];};
    struct shape_cmd {int type; int idx;};      // 0 box, 1 stl (in input order)
    std::vector<shape_box> boxes;
    std::vector<shape_stl> stls;
    std::vector<shape_cmd> shapes;
    std::vector<fix_box> fixes;
    std::vector<std::pair<std::string,Vec3>> mon_pts;

    std::vector<material> mats;

    // nodes
    std::vector<Vec3> X, x, v, fint, fext, fcon;
    std::vector<double> m, vnode;
    std::vector<unsigned char> fixed;           // bit 0,1,2: fixed dof
    std::vector<int> nalive;                    // number of intact elements per node
    std::vector<int> node_elem_start, node_elem; // node -> elements (CSR)
    std::vector<int> orphan;
    std::vector<double> Mcpl, mfl;
    std::vector<Vec3> Mu_cpl, fcpl;

    // elements
    std::vector<element> elems;
    std::vector<gpstate> gps;
    int ngp = 1;
    double dN0[8][3];                           // centre derivatives (reduced)
    double dNg[8][8][3];                        // [gp][node][dir] (full)
    double gam[4][8];                           // hourglass base vectors
    double Vel = 0.0, helem = 0.0;

    std::vector<face> faces;
    int surf_version = 0;
    bool surf_dirty = false;
    int nbodies = 1;

    std::vector<monitor> mons;
    coupling_options copt;

    // options
    int full_int = 0;
    double hg_coef = 0.1;
    double cfl = 0.5;
    double alpha_damp = 0.0;
    double relax_time = 0.0, relax_alpha = 0.0;
    double bulkq1 = 0.06, bulkq2 = 1.2;
    double erode_J = 0.05;
    bool plane_strain = false;
    bool ground_on = false;
    double zground = 0.0, kground = 1.0, mu_ground = 0.5;
    bool contact_on = true;
    double kcontact = 1.0, mu_contact = 0.5, contact_dist = 0.8, contact_zeta = 0.3;
    Vec3 grav = Vec3(0.0,0.0,-9.81);

    // state
    double t = 0.0;
    double dtcrit = 0.0;
    int nsub_last = 0;
    double wdiss = 0.0;
    bool built = false;
};

#endif
