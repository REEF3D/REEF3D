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

#ifndef FEM_SOLID_H_
#define FEM_SOLID_H_

// Explicit Total-Lagrangian finite element solver for solids
// (elastic deformation, stresses, material failure and collapse).
//
//  - mesh: 8-node hexahedra on a uniform voxel lattice, filled from boxes
//    and closed STL surfaces (shared lattice nodes glue touching bodies);
//    the surface nodes are snapped onto the STL / box surfaces, so only the
//    elements at the surface are distorted
//  - simple input: material presets (concrete C30, steel S355, ...), supports
//    in words (fix base / bed), element size from the fluid grid, gravity
//    settling, check mode with natural frequencies, engineering summary
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
#include<array>
#include<map>
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
        std::string name;               // preset name or type
        double rho = 1000.0, E = 1.0e6, nu = 0.3;
        double lambda = 0.0, mu = 0.0, cp = 0.0;
        // J2
        double sigy = 0.0, H = 0.0, epsfail = 0.0;
        // concrete
        double ft = 0.0, Gf = 0.0, fc = 0.0, Gc = 0.0;
        // erosion threshold for damage
        double derode = 0.99;
        // rigid: free bodies made only of rigid materials move as rigid bodies
        // (no stresses, no stiffness time step limit), e.g. floating debris
        bool rigid = false;
        bool Egiven = true;             // false: 'material rigid <rho>' (no elastic modulus given)
        // impact of rigid debris: effective stiffness [N/m] (0: axial bar E A / L)
        // and crushing force [N] that caps the contact force (0: no cap);
        // 2D: per metre width
        double kdebris = 0.0, fcrush = 0.0;
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
        int geo = -1;                   // geometry: -1 regular voxel, else index into geos
        int shape = -1;                 // input shape that created the element
        bool alive = true;
        bool rigid = false;             // part of a rigid body: no internal forces
        double svm = 0.0;               // von Mises (Cauchy) stress, averaged over the GPs
        double util = -1.0;             // utilisation (stress / strength), -1: elastic material
        double J = 1.0;
    };

    // reference geometry of a hex8 element
    struct egeom
    {
        double dN0[8][3];               // uniform (mean) gradient, one-point integration
        double gam[4][8];               // Flanagan-Belytschko hourglass vectors
        double bb;                      // sum of |dN0|^2
        double V;                       // volume
        double h;                       // crack-band length V^(1/3)
        double L;                       // characteristic length for the time step (V / max face area)
        double dNg[8][8][3];            // gradients at the 2x2x2 Gauss points
        double wg[8];                   // Gauss weights times det J
        double mw[8];                   // nodal mass weights (integral of N_a)
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
        int    loads = 0;               // 0: hybrid (parcels, corrected slowly toward the probed pressure),
                                        // 1: pressure integration (explicit), 2: attached fluid parcels only
        double hybrid_tau = 20.0;       // hybrid: filter time of the correction [fluid time steps]
        int    shear = 0;               // 1: tangential (no-slip) reaction of the parcels on the solid
        int    air_forcing = 0;         // 1: direct forcing also in the air (more than 1.5 cells from the water)
        int    walls = 1;               // 1: the boundaries of the fluid domain are walls for free bodies and debris
        double added_mass = 1.0;        // factor on the added-mass estimate of rigid bodies (stabilisation), 0: off
        int    settle = 1;              // 1: settle under gravity before the flow starts
        int    check = 0;               // 1: write the check report and stop
        double resolution = 1.0;        // element size in fluid cells if no lattice is given
        int    fix_bed = 0;             // 1: fix nodes touching the bed / solids of the fluid grid
    };

    // results of the check mode
    struct check_info
    {
        double mass = 0.0, volume = 0.0;
        Vec3 cog = Vec3::Zero(), lo = Vec3::Zero(), hi = Vec3::Zero();
        int nfixed = 0;
        double freq[3] = {0.0,0.0,0.0};     // Rayleigh estimate of the first frequency in x, y, z [Hz]
        bool freq_ok[3] = {false,false,false};
        double sw_maxdisp = 0.0, sw_maxvm = 0.0, sw_maxutil = -1.0;
        Vec3 sw_support = Vec3::Zero();
        bool settled = false;
        double thin_fraction = 0.0;         // fraction of surface elements with exposed opposite faces
        double surface_area = 0.0;
        std::vector<std::string> warnings;
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
    bool plane_strain_on() const {return plane_strain;}

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
    // structural damping: ratio of critical damping at the first natural
    // frequency, acts on the deformation only (rigid motion of free parts and
    // debris is not damped). prepare_damping() estimates the frequency.
    void set_damping_ratio(double z) {zeta = z;}
    void set_rigid_contact_speed(double c) {c_rigid = c;}
    void set_debris_damping(double z) {debris_zeta = z;}

    // contact of free bodies and debris with planes (domain walls) and with the
    // bed / solids of the fluid grid (signed distance phi and normal sampled per
    // node at the start of each fluid step, linearised during the substeps)
    void add_contact_plane(const Vec3& n,double d) {planes.push_back({n.normalized(),d});}
    void set_bed_sample(int i,double phi,const Vec3& n);
    void clear_bed_samples() {bed_ok.assign(nnode(),0); bed_phi.assign(nnode(),0.0); bed_n.assign(nnode(),Vec3::Zero()); bed_x.assign(nnode(),Vec3::Zero());}
    bool free_node(int i) const {return m[i]>0.0 && (body[i]<0 || !body_fixed[body[i]]);}
    bool bed_contact() const {return bed_on;}
    bool bed_sample(int i,double& phi,Vec3& n) const    // distance to the bed at the current position, bed normal
    {if(bed_ok.empty() || !bed_ok[i]) return false; phi = bed_phi[i] + bed_n[i].dot(x[i]-bed_x[i]); n = bed_n[i]; return true;}
    void set_bed_contact(bool b) {bed_on = b;}

    // rigid bodies: free bodies (no supports) made only of rigid materials
    struct rigid_body
    {
        std::vector<int> nodes;
        std::vector<Vec3> r0;           // node positions relative to the centre of mass, reference state
        double M = 0.0;                 // mass
        Vec3 c0 = Vec3::Zero(), c = Vec3::Zero(), V = Vec3::Zero(), L = Vec3::Zero(), w = Vec3::Zero();
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity(), I0 = Eigen::Matrix3d::Zero();
        double vmax = 0.0, dmax = 0.0, fcmax = 0.0;   // max speed, max displacement of the centre, max contact force
        double fcstep = 0.0;            // max contact force in the last fluid step
        double Vol = 0.0;               // volume
        double Aunit = 0.0, A = 0.0;    // added mass for the stabilisation: per unit fluid density, current
        Vec3 V0 = Vec3::Zero(), w0 = Vec3::Zero(), a_prev = Vec3::Zero(), al_prev = Vec3::Zero();
        // impact: effective contact stiffness [N/m] (debris acting as a spring), crushing force cap [N],
        // length of the axial bar, how k was obtained (0 given, 1 E A/L, 2 rigid_contact_speed)
        double k = 0.0, Fcap = 0.0, Lbar = 0.0;
        int ksrc = 0;
        Vec3 Jc = Vec3::Zero(), Hc = Vec3::Zero();   // contact impulse and its moment over the fluid step
    };
    void set_rigid_added_mass(int k,double A) {rbs[k].A = A;}
    int rigid_of_node(int i) const {return rnode.empty() ? -1 : rnode[i];}
    int n_rigid() const {return (int)rbs.size();}
    const rigid_body& rigid(int k) const {return rbs[k];}
    bool is_rigid_node(int i) const {return !rnode.empty() && rnode[i]>=0;}
    bool node_rigid_material(int i) const       // all elements at the node are of rigid materials
    {for(int q=node_elem_start[i]; q<node_elem_start[i+1]; ++q) if(!mats[elems[node_elem[q]].mat].rigid) return false; return node_elem_start[i]<node_elem_start[i+1];}
    int n_deformable() const {int k = 0; for(const element& e : elems) if(e.alive && !e.rigid) ++k; return k;}
    bool rigid_material_deformable() const      // rigid material in a supported or mixed body
    {for(const element& e : elems) if(e.alive && mats[e.mat].rigid && !e.rigid) return true; return false;}
    double damping_ratio() const {return zeta;}
    void prepare_damping();
    double damping_alpha() const {return alpha_struct;}
    double damping_frequency() const {return f_damp;}
    // first natural frequency per direction (Rayleigh quotient of the static
    // deflection under 1 g), the state is restored afterwards
    void rayleigh_frequencies(double f[3],bool ok[3]);
    void set_ground(double z,double kfac,double mu) {ground_on = true; zground = z; kground = kfac; mu_ground = mu;}
    void set_contact(bool on,double kfac,double mu) {contact_on = on; kcontact = kfac; mu_contact = mu;}
    void set_snap(bool on) {snap_on = on;}
    bool lattice_given() const {return hx>0.0 && hy>0.0 && hz>0.0;}
    void set_default_spacing(double h,double hy=-1.0);   // lattice spacing when the input has none (resolution)
    bool ground() const {return ground_on;}
    static bool preset(const std::string& type,const std::string& name,material&);   // material presets
    static std::string preset_list();

    // ------------------------------------------------------------------
    // tools
    // ------------------------------------------------------------------
    // static equilibrium under gravity (kinetic damping), velocities zero afterwards
    bool settle(int maxsteps=200000,double tol=1.0e-4,double* residual=nullptr,bool supported_only=false,bool elastic=false);
    // check: geometry, mass, supports, frequencies, self-weight response, warnings
    check_info check();
    void write_check(std::ostream&,const check_info&) const;
    // fix the nodes of the given list (all dofs)
    void fix_nodes(const std::vector<int>& nodes);
    // utilisation (stress / strength) of every intact element from the current state
    void update_utilisation();
    double max_utilisation(int* elem=nullptr) const;
    double max_damage(int* elem=nullptr) const;
    double max_plastic_strain(int* elem=nullptr) const;
    Vec3 elem_centre(int e) const;
    Vec3 support_moment() const;                // overturning moment of the support forces about the base centre
    Vec3 base_centre() const {return base_c;}
    double eroded_mass_fraction() const;
    bool has_supports() const;
    int  material_count() const {return (int)mats.size();}
    const material& mat(int i) const {return mats[i];}
    double min_density() const;

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
    const Vec3& vel_mean(int i) const {return vbar[i];}     // mean velocity over the last step (coupling)
    double mass(int i) const {return m[i];}
    double node_volume(int i) const {return vnode[i];}
    double node_density(int i) const {return vnode[i]>0.0 ? m[i]/vnode[i] : 0.0;}
    void set_vel(int i,const Vec3& vv) {v[i] = vv;}
    bool is_fixed(int i) const {return fixed[i]!=0;}
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
    void add_coupling(int i,const Vec3& M,const Vec3& Mu) {Mcpl[i] += M; Mu_cpl[i] += Mu;}
    void set_fluid_mass(int i,double mf) {mfl[i] = mf;}
    // fluid load of the parcels and the enclosed fluid in the last step
    Vec3 parcel_load(int i) const {return fcpl[i] - mfl[i]*grav;}
    // node belongs to an intact body with supports (not a free fragment, not debris)
    bool on_supported_body(int i) const {return body[i]>=0 && body_fixed[body[i]];}
    // weight of the final RK stage of the flow solver: the forcing in that stage
    // only removes alpha times the momentum the fluid exchanges over a step
    void set_coupling_alpha(double a) {alpha_cpl = a;}
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

    void snap_surface();
    void make_fixes();
    bool element_geometry(const Vec3* Xa,egeom& G) const;   // false if inverted
    const egeom& geom(const element& el) const {return el.geo<0 ? vgeo : geos[el.geo];}

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
    struct shape_stl {std::string file; int mat; bool remove; double scale = 1.0, rot = 0.0; Vec3 move = Vec3::Zero();};
    struct fix_box {double x0,x1,y0,y1,z0,z1; bool f[3]; int mode = 0;};   // mode 0 box, 1 base, 2 top
    struct mon_auto {std::string name; int mode;};
    struct shape_cmd {int type; int idx;};      // 0 box, 1 stl (in input order)
    std::vector<shape_box> boxes;
    std::vector<shape_stl> stls;
    std::vector<shape_cmd> shapes;
    std::vector<fix_box> fixes;
    std::vector<std::pair<std::string,Vec3>> mon_pts;
    std::vector<mon_auto> mon_autos;
    std::vector<int> vox_shape;                 // voxel -> shape that filled / carved it (-1 none)
    std::vector<std::vector<std::array<double,9>>> stl_tris;   // triangles per STL shape (transformed)
    int cur_mat = -1;                           // material id used by shapes without id

    std::vector<material> mats;

    // nodes
    std::vector<Vec3> X, x, v, fint, fext, fcon;
    std::vector<double> m, vnode;
    std::vector<unsigned char> fixed;           // bit 0,1,2: fixed dof
    std::vector<int> nalive;                    // number of intact elements per node
    std::vector<int> node_elem_start, node_elem; // node -> elements (CSR)
    std::vector<int> orphan;
    std::vector<double> mfl;
    std::vector<Vec3> Mcpl;                     // attached fluid mass per direction
    std::vector<Vec3> Mu_cpl, fcpl, vbar;

    // elements
    std::vector<element> elems;
    std::vector<gpstate> gps;
    int ngp = 1;
    egeom vgeo;                                 // regular voxel
    std::vector<egeom> geos;                    // distorted (snapped) elements
    double Vel = 0.0, helem = 0.0;
    Vec3 base_c = Vec3::Zero();

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
    double zeta = 0.02, alpha_struct = 0.0, f_damp = 0.0;
    double c_rigid = 40.0;                      // rigid materials without E: contact stiffness from E = rho c^2
    struct cplane {Vec3 n; double d;};          // n.x >= d is outside the wall
    std::vector<cplane> planes;
    bool bed_on = false;
    std::vector<unsigned char> bed_ok;
    std::vector<double> bed_phi;
    std::vector<Vec3> bed_n, bed_x;
    void contact_surface(int i,double pen,const Vec3& n,double mu);
    // contact of rigid bodies: every rigid body acts as one spring of its effective
    // stiffness k against each partner (wall, bed, another body, the deformable parts)
    struct rpair {int a, b; long long key; Vec3 n; double pen, vn;};
    std::vector<rpair> rpairs;
    std::map<std::pair<long long,long long>,double> rset;   // permanent set of crushed contact points (group, point)
    void rigid_contact_forces();
    double debris_zeta = 0.5;                   // damping ratio of the debris contact (unloading only)
    std::vector<Vec3> tspring;                  // tangential contact spring of the nodes on walls / bed / ground
    std::vector<unsigned char> touched;
    double dts_cur = 0.0;
    double cp_contact = 0.0;                    // wave speed of the contact penalty (largest)
    std::vector<double> cnode;                  // contact wave speed per node
    std::vector<rigid_body> rbs;
    std::vector<int> rnode;                     // node -> rigid body (-1: deformable)
    void setup_rigid();
    bool bodies_near() const;
    void rigid_step(double dts,const std::vector<Vec3>& F);
    std::vector<int> body;                      // node -> connected body (-1: debris / no intact element)
    std::vector<unsigned char> body_fixed;      // body has supported nodes
    void rigid_velocity(std::vector<Vec3>& vr) const;   // rigid-body velocity of the free bodies
    double relax_time = 0.0, relax_alpha = 0.0;
    double bulkq1 = 0.06, bulkq2 = 1.2;
    double erode_J = 0.05;
    bool plane_strain = false;
    bool ground_on = false;
    double zground = 0.0, kground = 1.0, mu_ground = 0.5;
    bool contact_on = true;
    bool snap_on = true;
    double alpha_cpl = 1.0;
    double kcontact = 1.0, mu_contact = 0.5, contact_dist = 0.8, contact_zeta = 0.3;
    Vec3 grav = Vec3(0.0,0.0,-9.81);

    // state
    double t = 0.0;
    double dtcrit = 0.0, dtcrit_el = 1.0e30, dtcrit_rig = 1.0e30;   // overall, elements, rigid-body impact
    bool rigid_contact_near(double dt) const;
    int nsub_last = 0;
    double wdiss = 0.0;
    bool built = false;
};

#endif
