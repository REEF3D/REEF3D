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

// Impermeable membrane (closed flexible fish cage) for REEF3D::NHFLOW.   ctrl.txt: X 330 1, A 520 1
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
//    correction as the mobility beta = 1/(1 + a K_n H) (d->MBETA; nhflow_poisson, nhflow_pjm).
//    Without it the correction -a/rho grad(q) puts O(a dq/(rho delta)) normal velocity through the
//    membrane in every stage. Only the full projection A 520 1 is supported.
//
// 3. Free surface. The depth-jump dissipation of the HLL/HLLC continuity flux is scaled with beta,
//    otherwise it moves mass through the layer although the velocity there vanishes.
//
// 4. Bag floor. NHFLOW has one free surface per column, so the columns below the floor see the
//    inner level in their hydrostatic pressure. The static overpressure -rho g dh below the floor
//    is prescribed (static_pressure_nhflow) and the continuity dissipation is switched off where it
//    acts (d->MCHI), so the solved non-hydrostatic pressure stays smooth.
//
// 5. Loads. Evaluated after the projection from the final velocity, F = rho H A (u^{n+1} - u_m) dV,
//    per triangle. In steady state they integrate to the pressure jump times the panel area.
//
// This first version carries prescribed (fixed) membranes; the node velocities xdot_ enter the
// forcing already, the structural solver is connected in a follow-up.

#include"net.h"
#include"increment.h"
#include<vector>
#include<array>
#include<string>
#include<Eigen/Dense>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class slice;

using namespace std;

struct membrane_param
{
    string name;
    int shape=0;                // 1: box, 2: cylinder
    double x0=0.0,x1=0.0,y0=0.0,y1=0.0;     // box footprint
    double xc=0.0,yc=0.0,R=0.0;             // cylinder footprint
    double zb=0.0,zt=0.0;                   // bottom (floor) and top of the bag
    double Rn=1.0e4, Rt=0.0;                // hydraulic resistance normal / tangential [m/s]
    double delta=-1.0;                      // half width of the smeared layer [m], <0: 1.5 max(dx,dy,dz)
    double h=-1.0;                          // target triangle edge length [m], <0: min cell size
    double fill=0.0;                        // initial inner water level above the undisturbed level [m]
    double printdt=-1.0;                    // vtp print interval [s], <0: no vtp output
    int poisson=1;                          // 1: membrane mobility in the pressure Poisson equation, 0: off (diagnostics)
};

class net_membrane final : public net, public increment
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

    // NHFLOW coupling, called per RK stage: mobility (builds the cell map), forcing, reaction
    void mobility_nhflow(lexer*, fdm_nhf*, ghostcell*, double);
    void forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void static_pressure_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void reaction_nhflow(lexer*, fdm_nhf*, ghostcell*, double, slice&, bool);

    // initial inner water level
    void fill_nhflow(lexer*, fdm_nhf*, ghostcell*);

    const membrane_param& param() const {return prm;}

private:
    // geometry
    void mesh(lexer*);
    void add_panel(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, int, int, const Eigen::Vector3d&, int);
    void add_cylinder_wall(int, int);
    void add_disk(int, int);
    void add_tri(int, int, int, const Eigen::Vector3d&, int);
    bool inside_footprint(double, double, double) const;
    bool outside_footprint(double, double, double) const;

    // cell map
    void build_map(lexer*, fdm_nhf*, ghostcell*);
    static Eigen::Vector3d closest_point(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&,
                                         const Eigen::Vector3d&, double&, double&, double&);
    Eigen::Vector3d membrane_vel(int, double, double, double) const;
    double indicator(double) const;
    double smoothstep(double) const;
    double footprint_distance(double, double, double&, double&) const;

    // output
    void print_timeseries(lexer*, fdm_nhf*, ghostcell*);
    void print_vtp(lexer*);

    int nMem;
    membrane_param prm;
    EigenMat empty_;

    // triangulation
    vector<Eigen::Vector3d> x_, xdot_;
    vector<array<int,3> > tri_;
    vector<Eigen::Vector3d> tn_, tc_;       // outward unit normal, centroid
    vector<double> ta_;                     // area
    vector<int> ttag_;                      // 0: wall, 1: floor
    double Afloor;

    // smeared layer
    double delta, Kn, Kt;

    struct cellentry
    {
        int i,j,k;
        int t1,t2;
        double d1,d2;
        double H1,H2;
        double w0,w1,w2;    // barycentric weights of the closest point on t1
    };
    vector<cellentry> cells_;
    vector<int> slot_;
    vector<double> xc_, yc_;                // local cell centres for the bounding box search

    // loads
    vector<double> tf_;                     // 3 per triangle, reaction on the membrane [N]
    double Fx,Fy,Fz,Fzfloor,Qleak,urelmax;
    double dh;                              // eta_in - eta_out used for the static floor pressure

    // print
    double printtime;
    int printcount;
    string outdir;
};

#endif
