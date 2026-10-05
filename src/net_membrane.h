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

#ifndef NET_MEMBRANE_H_
#define NET_MEMBRANE_H_

// Impermeable membrane (closed flexible fish cage) for REEF3D::NHFLOW.   ctrl.txt: X 330 1, A 520 1 or 2
// Geometry and parameters in membrane.dat (see net_interface_membrane_nhflow.cpp).
//
// 1. Porous jump. The membrane is a triangulated surface. On the fluid side it is a thin, anisotropic
//    porous layer of half width delta around the surface:
//
//      f = -H(d) A (u - u_m),     A = K_n (n n^T) + K_t (I - n n^T),     K = R/(1.5 delta)
//
//    H(d) = 1 for d < delta/2 with a cosine taper to 0 at d = delta; its integral across the layer
//    is 1.5 delta, so the pressure jump over the layer is  dp/rho = R_n u_n  (hydraulic resistance R
//    in m/s, leakage velocity u_n = (dp/rho)/R_n). The plateau spreads the pressure / free surface
//    drop linearly over the layer instead of concentrating it in one cell. The resistance is
//    integrated implicitly over the stage, (u - u_m) <- (I + a A)^{-1} (u - u_m), a = alpha dt.
//    Cells near an edge or corner (floor - wall) get the normal resistance of the closest triangle
//    and of the closest triangle with a different normal, so no leakage path opens around corners.
//
// 2. Projection. The same implicit factor enters the pressure Poisson equation and the velocity
//    correction as the mobility beta = 1/(1 + a K_n H) (d->MBETA; nhflow_poisson(_pcorr),
//    nhflow_pjm(_corr), nhflow_membrane_beta.h): in the matrix per face, in the collocated correction
//    as the distance-weighted mean of the two compact face gradients with their face mobilities.
//    A 520 2 (incremental): the old pressure gradient of the predictor is taken out before the forcing
//    and put back through the same operator (net_interface::membrane_pgrad), so P^n and the increment
//    act like one full pressure. Optional further projection passes per stage (projections n).
//
// 3. Free surface. At faces next to the membrane or below the bag floor the continuity flux uses the
//    Rhie-Chow face velocity of the last projection (with projections > 1 the central face velocity of
//    the converged wide divergence), so the free surface follows the projected velocity field and no
//    mass moves through the membrane layer or past the floor corner that the projection does not see.
//
// 4. Bag floor. NHFLOW has one free surface per column, so the columns below the floor see the
//    inner level in their hydrostatic pressure. The static overpressure below the floor is prescribed
//    (static_pressure_nhflow, floorpressure 0/1/3), so the solved non-hydrostatic pressure stays smooth.
//    For a moving membrane the footprint follows the floor edge and the floor height is taken per
//    column from the floor triangles.
//
// 5. Loads. Evaluated after the projection from the final velocity, F = rho H A (u^{n+1} - u_m) dV,
//    per triangle and, with the barycentric weights of the closest point, per node. In steady state
//    they integrate to the pressure jump times the panel area.
//
// 6. Structure (net_membrane_structure.cpp), membrane.dat 'structure':
//      fixed      the membrane stays in place (default)
//      rigid      the whole membrane moves with the floating body (X 10), loads go to the body
//      flexible   mass-spring membrane (edge springs from E t, edge dampers, fabric weight and
//                 buoyancy, sinker weight along the floor edge). Like the nets (X 320), the nodes at
//                 the top edge ('attach') follow the floating body, or stay in place without X 10, and
//                 the force of the membrane on them goes to the body as an external force.
//    The flexible membrane is advanced once per time step in the final RK stage with a linearised backward
//    Euler step, in which the porous-jump load is taken implicitly in the node velocity (the fluid of the
//    layer is tied to the membrane like a stiff damper, which a lagged, explicit coupling can not carry).
//    Stable, equilibrium unchanged, but the fast dynamics are damped (see net_membrane_structure.cpp).
//    The fluid layer of a moving membrane moves with it (R_t = R_n by default): the Poisson mobility is
//    isotropic, so the layer fluid can not be moved along the membrane by pressure.
//
// 7. Strong coupling (membrane.dat 'coupling iterated', flexible membranes; net_membrane_coupling.cpp). The
//    projection of every RK stage is repeated until the node velocities the fluid used and the ones the structure
//    returns for the resulting loads agree (nhflow_forcing::projection). The structure takes a backward-Euler step
//    of the stage with its real inertia; the local response of the layer fluid serves as Robin preconditioner
//    and drops out at convergence; IQN-ILS accelerates the fixed point (iqn_ils.h). No artificial inertia.
//    Defaults: membrane mesh 1.5 delta (shorter structural modes are not seen by the layer and stall the iteration),
//    Robin matrix x16 with the tolerances / 16.
//
// 8. Current. In a steady current at low Froude number use the HLLC flux (A 511 2): HLL dissipates the shear at the
//    edge of the frozen membrane layer with the gravity wave speed, and the thick, friction-like wake gave drag
//    coefficients 4 - 7 times too high for the cage of Strand et al. (2013) (U/sqrt(g h) ~ 0.03); HLLC carries the
//    shear wave. floorpressure 3 is unstable in a current, use 0. Slack (partly filled) bags: coupling staggered with
//    a lower resistance (1e3) and compression 0.001 (the iterated coupling does not converge for wrinkling fabric).
//    Shapes: box, cylinder, cylcone (cylinder on a cone; the cone is the floor, sloped floor in the static pressure).

#include"net.h"
#include"increment.h"
#include"vtp3D.h"
#include<vector>
#include<array>
#include<string>
#include<Eigen/Dense>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class slice;

#include<functional>
using namespace std;

struct membrane_param
{
    string name;
    int shape=0;                // 1: box, 2: cylinder, 3: cylinder on a cone (cylcone)
    double x0=0.0,x1=0.0,y0=0.0,y1=0.0;     // box footprint
    double xc=0.0,yc=0.0,R=0.0;             // cylinder footprint
    double zb=0.0,zt=0.0;                   // bottom (floor) and top of the bag; cylcone: zb the cone tip
    double zc=0.0;                          // cylcone: cone base = bottom of the cylinder
    double Rn=1.0e4, Rt=-1.0;               // hydraulic resistance normal / tangential [m/s]; Rt<0: 0 fixed, Rn moving
                                            // (link mode: 0)
    double delta=-1.0;                      // half width of the smeared layer [m], <0: 1.5 max(dx,dy,dz)
    double h=-1.0;                          // target triangle edge length [m], <0: min cell size (fixed), max cell size
                                            // (moving), 1.5 delta (coupling iterated)
    double fill=0.0;                        // initial inner water level above the undisturbed level [m]
    double filling=-1.0;                    // filling level V_water/V_bag (<0: not given, fill applies)
    double drain=0.0;                       // >0: the missing water is pumped out of the bag over this time [s]
    double printdt=-1.0;                    // vtp print interval [s]; <0: NHFLOW print control (P 20 / P 30), 0: off
    int projections=1;                      // projection passes per stage (1: Rhie-Chow flux, >1: converged wide divergence)
    int poisson=1;                          // 1: membrane mobility in the pressure Poisson equation, 0: off (diagnostics)
    int link=0;                             // pressure coupling: 0 layer (isotropic mobility of the layer cells), 1 link
                                            // (porous-jump mobility only on the links crossing the membrane)
    int floorp=-1;                          // static pressure below the floor: 0 uniform dh, 1 local, 3 averaged ramp;
                                            // -1: 3 for a fixed membrane, 0 for a moving one
    double tau=2.0;                         // averaging time of the ramp shape (floorpressure 3), frozen at 2 tau [s]

    // structure
    int structure=0;                        // 0 fixed, 1 rigid (moves with the floating body), 2 flexible
    double mA=1.0;                          // fabric mass per area [kg/m^2]
    double rhom=1300.0;                     // fabric density [kg/m^3] (buoyancy)
    double EA=5.0e5;                        // membrane stiffness E t [N/m]
    double zeta=0.1;                        // damping ratio of the edge dampers
    double compr=0.01;                      // stiffness in compression as a fraction of E t (wrinkling fabric)
    double sinker=0.0;                      // submerged weight along the floor edge [N/m]
    double zattach=-1.0e20;                 // nodes at or above this height are attached, default: top edge
    double Mbody=-1.0;                      // added mass of the coupling to the floating body [kg], <0: 2 rho V_bag

    // fluid-structure coupling of a flexible membrane
    int coupling=-1;                        // 0 staggered (once per step, implicit porous damper), 1 iterated per stage,
                                            // -1 not given: iterated for a flexible membrane
    double crtol=1.0e-3;                    // iterated: relative tolerance of the node velocities
    double catol=1.0e-5;                    //           absolute tolerance [m/s] (rms)
    int citer=50;                           //           maximum iterations per stage
    int creuse=8;                           //           IQN-ILS: converged stages kept for the next ones
    double crelax=0.5;                      //           relaxation of the first iteration without history
    int clog=0;                             //           1: residual of every iteration on screen
    double crobin=16.0;                     //           scale of the Robin matrix (local layer response)
    int ccols=100;                          //           IQN-ILS: maximum number of columns
    double cfilt=1.0e-2;                    //           IQN-ILS: QR filter (drop columns nearly dependent on newer ones)
    int cqn=0;                              //           quasi-Newton: 0 IQN-ILS, 1 IQN-IMVJ

    // flexible collar (membrane.dat 'collar'): the top edge of a flexible bag is a floating pipe ring
    int collar=0;
    double cD=0.5;                          // pipe diameter [m]
    double cm=60.0;                         // mass per length [kg/m] (pipe, brackets, walkway)
    double cEA=4.5e7;                       // axial stiffness [N]
    double cEI=1.2e6;                       // bending stiffness [N m^2]
    double cCd=1.0, cCa=1.0;                // Morison drag and added-mass coefficients (normal to the pipe axis)
    vector<array<double,5> > moor;          // mooring springs: anchor x y z, stiffness [N/m], pretension [N]
};

class net_membrane final : public net, public increment, private vtp3D
{
public:
    net_membrane(int, const membrane_param&);
    virtual ~net_membrane();

    // net interface (the membrane runs through the membrane_* calls of net_interface)
    void start_cfd(lexer*, fdm*, ghostcell*, double, Eigen::Matrix3d&, bool) override final {};
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, double, Eigen::Matrix3d&, bool) override final {};
    void initialize_cfd(lexer*, fdm*, ghostcell*) override final {};
    void initialize_nhflow(lexer*, fdm_nhf*, ghostcell*) override final;
    void netForces(lexer*, double&, double&, double&, double&, double&, double&) override final;

    const EigenMat& getLagrangePoints() override final {return empty_;}
    const EigenMat& getLagrangeForces() override final {return empty_;}
    const EigenMat& getCollarVel() override final {return empty_;}
    const EigenMat& getCollarPoints() override final {return empty_;}

    // NHFLOW coupling, called per RK stage: kinematics, mobility (builds the cell map), forcing, reaction
    void kinematics_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void mobility_nhflow(lexer*, fdm_nhf*, ghostcell*, double);
    void forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void static_pressure_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void reaction_nhflow(lexer*, fdm_nhf*, ghostcell*, double, slice&, bool);

    // floating body: body frame at the start, current kinematics, load on the body
    void attach_body(lexer*, const Eigen::Vector3d&, const Eigen::Matrix3d&);
    void set_body(const Eigen::Vector3d&, const Eigen::Matrix3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
    void body_load(lexer*, const Eigen::Vector3d&, const Eigen::Matrix3d&, double&, double&, double&, double&, double&, double&) const;
    double body_addedmass(lexer*) const;
    Eigen::Matrix3d body_addedinertia(lexer*) const;
    Eigen::Matrix3d body_stiffness() const {return (moving() && body_) ? Eigen::Matrix3d(Jb_.block<3,3>(0,0)) : Eigen::Matrix3d::Zero();}

    // initial inner water level, or water pumped out of the bag over the drain time
    void fill_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void drain_nhflow(lexer*, fdm_nhf*, ghostcell*);

    // strong coupling (net_membrane_coupling.cpp): forcing update for new node velocities, coupling iteration
    bool iterated() const {return prm.coupling==1 && prm.structure==2;}
    void reforce_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    bool couple_nhflow(lexer*, fdm_nhf*, ghostcell*, int, double, slice&, int);

    const membrane_param& param() const {return prm;}

private:
    // geometry
    void mesh(lexer*);
    void merge_nodes();
    void update_geometry();
    void add_panel(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, int, int, const Eigen::Vector3d&, int);
    void add_cylinder_wall(int, int, double);
    void add_disk(int, int, double, double);
    void add_tri(int, int, int, const Eigen::Vector3d&, int);
    bool inside_footprint(double, double, double) const;
    bool outside_footprint(double, double, double) const;

    // cell map
    void build_map(lexer*, fdm_nhf*, ghostcell*);
    static Eigen::Vector3d closest_point(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&,
                                         const Eigen::Vector3d&, double&, double&, double&);
    Eigen::Vector3d membrane_vel(int, double, double, double) const;
    Eigen::Vector3d membrane_vel(const vector<Eigen::Vector3d>&, int, const double*) const;
    void compute_loads(lexer*, fdm_nhf*, ghostcell*, slice&);
    double indicator(double) const;
    double dindicator(double) const;
    double smoothstep(double) const;
    double footprint_distance(double, double, double&, double&) const;
    void footprint_weight(double, double, double&, double&, double&) const;
    void footprint_weight_ext(double, double, double&, double&, double&) const;
    void footprint_weight_int(double, double, double&, double&, double&) const;
    double dstep(double) const;
    void floor_geometry(lexer*);
    double zfloor(lexer*, int, int) const;
    double zbag() const {return prm.shape==3 ? prm.zc - (prm.zc-prm.zb)/3.0 : prm.zb;}   // flat floor of the same volume
    double Abag(lexer*) const;              // waterplane area of the bag
    bool sloped() const {return prm.shape==3;}  // floor height varies over the footprint

    // structure (net_membrane_structure.cpp)
    void ini_structure(lexer*, ghostcell*);
    void advance_structure(lexer*, double);
    void structure_solve(lexer*, double, int, const vector<Eigen::Matrix3d>&, const vector<Eigen::Vector3d>&,
                         const vector<Eigen::Vector3d>&, const vector<Eigen::Vector3d>&, bool=false);
    void free_factor();
    struct structure_factor *sfac_=nullptr;   // factorised structure matrix, reused in the coupling iterations of a stage
    bool sfac_ok_=false;
    void internal_forces(lexer*, const vector<Eigen::Vector3d>&, const vector<Eigen::Vector3d>&, vector<Eigen::Vector3d>&) const;
    void external_forces(lexer*, vector<Eigen::Vector3d>&) const;
    void update_body_load(lexer*);
    void attach_response(lexer*, double, const vector<Eigen::Matrix3d>* =nullptr);
    void update_body_fluid_load(lexer*);
    Eigen::Vector3d attached_position(int) const;
    Eigen::Vector3d body_velocity(const Eigen::Vector3d&) const;
    bool moving() const {return prm.structure>0;}

    // flexible collar (net_membrane_collar.cpp)
    void ini_collar(lexer*);
    void sample_collar(lexer*, fdm_nhf*, ghostcell*);
    void collar_step_end(lexer*);
    void collar_assemble(lexer*, double, const vector<int>&, const function<void(int,int,const Eigen::Matrix3d&)>&, Eigen::VectorXd&);
    void collar_output(lexer*);
    bool collar() const {return prm.collar==1 && prm.structure==2;}
    vector<int> cn_;                        // collar nodes in ring order
    vector<double> cl_, cL0_;               // length per node, rest length of the segment i -> i+1
    vector<Eigen::Vector3d> cxr_;           // reference ring relative to its centroid
    vector<Eigen::Vector3d> cuf_, cuf0_, caf_;  // fluid velocity at the collar nodes (stage, last step), fluid acceleration
    vector<double> ceta_;                   // free surface (absolute) above the collar nodes
    double ctn_=-1.0;                       // time of cuf0_
    vector<int> mfl_;                       // fairlead node of each mooring spring
    vector<double> mL0_;                    // initial length of each mooring spring
    vector<double> mT_;                     // mooring tensions
    double cMmax_=0.0, cNmax_=0.0;          // largest bending moment and axial force in the collar

    // output (net_membrane_print.cpp)
    void print_timeseries(lexer*, fdm_nhf*, ghostcell*);
    void print_vtp(lexer*);
    bool print_now(lexer*);

    int nMem;
    membrane_param prm;
    EigenMat empty_;

    // triangulation
    vector<Eigen::Vector3d> x_, xdot_;
    vector<array<int,3> > tri_;
    vector<Eigen::Vector3d> tn_, tc_;       // outward unit normal, centroid
    vector<Eigen::Vector3d> tn0_;           // outward unit normal of the undeformed mesh (orientation classes)
    vector<double> ta_;                     // area
    vector<int> ttag_;                      // 0: wall, 1: floor
    double Afloor;

    // smeared layer
    double delta, Kn, Kt;
    double rmap_=-1.0;                      // radius of the cell map (delta; link mode: at least the longest link)

    struct cellentry
    {
        int i,j,k;
        int ns;             // distinct surface orientations within delta (up to 3)
        int t[3];           // closest triangle of each orientation
        double dd[3], H[3];
        double bw[3][3];    // barycentric weights of the closest point on t[r]
        double Hmax;
        int tc;             // closest triangle overall
        double dc, Hc;
        double w0,w1,w2;    // barycentric weights of the closest point on tc
    };
    vector<cellentry> cells_;
    vector<int> slot_;
    vector<double> xc_, yc_;                // local cell centres for the bounding box search

    // link mode (net_membrane_link.cpp): cell centres and nodes on either side of the membrane, link mobilities,
    // loads from the forcing impulse and the pressure difference across the blocked links
    void link_ini(lexer*);
    void pseudonormals();
    int side_of(const Eigen::Vector3d&, int, double, double, double) const;
    int side_near(const Eigen::Vector3d&, const cellentry&) const;
    void link_mobility(lexer*, fdm_nhf*, ghostcell*, double);
    void link_loads(lexer*, fdm_nhf*, ghostcell*, slice&, const function<void(int,const double*,const Eigen::Vector3d&)>&);
    double *sideC_=nullptr;                 // side of the cell centres in the layer: +1 outside, -1 inside, 0 not set
    vector<signed char> sideN_;             // side of the nodes (FIJK) in the layer
    vector<int> sideCq_, sideNq_;           // indices set in this stage
    vector<Eigen::Vector3d> vpn_, epn_;     // angle-weighted vertex and edge pseudo normals
    vector<array<int,2> > etri_;            // triangles of each edge (-1: boundary edge)
    vector<Eigen::Vector3d> fimp_;          // per layer cell: force of the forcing on the membrane in this stage [N]
    struct blockedlink {int i,j,k,dir,t; double w[3];};
    vector<blockedlink> blocked_;           // links crossing the membrane in this stage


    // loads
    vector<double> tf_;                     // 3 per triangle, reaction on the membrane [N]
    vector<double> nf_;                     // 3 per node, reaction on the membrane [N]
    double Fx,Fy,Fz,Fzfloor,Qleak,urelmax;
    double dh, etaref;                      // eta_in - eta_out and eta_out used for the static floor pressure
    vector<double> etab_;                   // averaged free surface (floorpressure 3)
    double erefb=0.0, dhb=0.0;
    double dVdrain_=0.0, Vdrained_=0.0;     // volume to pump out of the bag, pumped so far

    // moving floor: footprint offset (floor edge centroid), floor height per column
    vector<int> ring_;                      // nodes on the floor edge
    Eigen::Vector3d ringc0_;
    double offx_=0.0, offy_=0.0, zring_=0.0;
    vector<double> zf_, zfx_, zfy_;         // floor height per column and its slope

    // structure
    vector<Eigen::Vector3d> x0_;            // initial node positions
    vector<Eigen::Vector3d> vs_;            // structural node velocity
    vector<Eigen::Vector3d> xan_;           // attached node positions at the start of the step
    vector<Eigen::Vector3d> rb_;            // node positions in the body frame
    vector<Eigen::Vector3d> nn_;            // node normals
    vector<array<int,2> > edge_;
    vector<double> L0_, ke_, ce_;
    vector<double> an_, mn_, ws_;           // node area, mass, sinker weight
    vector<char> att_;                      // attached to the floating body (or fixed)
    vector<int> grp_;                       // representative node of the group (2D: the node pairs across the slice)
    vector<vector<int> > gm_;               // groups of free nodes, moving together
    vector<array<int,3> > tedge_;           // edges of each triangle
    vector<int> flo_;                       // floor nodes
    vector<char> cornerE_, cornerN_;        // edges between panels of different orientation (floor edge, box corners), their nodes
    int nsub_=0;
    double Tmax_=0.0, vmax_=0.0;

    // floating body
    bool body_=false;
    Eigen::Vector3d cb_=Eigen::Vector3d::Zero(), vb_=Eigen::Vector3d::Zero(), wb_=Eigen::Vector3d::Zero();
    Eigen::Matrix3d Rb_=Eigen::Matrix3d::Identity();
    Eigen::Vector3d Fb_=Eigen::Vector3d::Zero(), Mb_=Eigen::Vector3d::Zero();   // load on the body, moment about the origin
    Eigen::Matrix3d Iab_=Eigen::Matrix3d::Zero();                           // 2 rho x inertia of the bag water about the CoG, body frame
    Eigen::Matrix<double,6,6> Jb_=Eigen::Matrix<double,6,6>::Zero();            // d(F, M_origin)/d(body translation, rotation)
    Eigen::Vector3d Ffl_=Eigen::Vector3d::Zero(), Mfl_=Eigen::Vector3d::Zero(); // fluid load on the attached edge (flexible), every stage
    Eigen::Vector3d cbn_=Eigen::Vector3d::Zero();                               // body position and orientation of Fb_, Mb_
    Eigen::Matrix3d Rbn_=Eigen::Matrix3d::Identity();

    // strong coupling (net_membrane_coupling.cpp)
    void layer_matrix(const cellentry&, const vector<Eigen::Vector3d>&, Eigen::Matrix3d&, Eigen::Vector3d&) const;
    void robin_matrix(lexer*, ghostcell*, double, slice&);
    void stage_begin(lexer*, ghostcell*, int, double, slice&);
    void stage_commit(lexer*);
    void group_vector(const vector<Eigen::Vector3d>&, Eigen::VectorXd&, bool) const;
    void group_scatter(const Eigen::VectorXd&, vector<Eigen::Vector3d>&) const;
    vector<Eigen::Vector3d> xn_, vsn_;      // structure at the start of the time step
    vector<Eigen::Vector3d> xk_, vsk_;      // structure after the last converged stage
    vector<Eigen::Vector3d> xbase_, vbase_; // base state of the current stage, (1-alpha) n + alpha k
    vector<Eigen::Vector3d> xdotf_;         // node velocities of the forcing of this stage
    vector<Eigen::Vector3d> xr_, vr_;       // structure result of the current iteration
    vector<Eigen::Matrix3d> Cr_;            // Robin matrix per node: local response of the layer fluid
    class iqn_ils *pqn_=nullptr;
    bool pending_=false;                    // converged stage result waiting for reaction_nhflow
    int citstep_=0, citmax_=0, cwarn_=0;
    double cres_=0.0, rbest_=0.0;

    // print
    double printtime;
    int printcount;
    string outdir, vtpdir;
    vector<double> pvdtime_;
};

#endif
