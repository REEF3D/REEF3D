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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#ifndef FSI_STRIP_H_
#define FSI_STRIP_H_

#include<vector>
#include <Eigen/Dense>
#include"beam.h"
#include"field4.h"

class lexer;
class fdm;
class ghostcell;
class field;
class turbulence;

using namespace std;

class fsi_strip : public beam, public increment
{
public:
	
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    typedef Eigen::Matrix<double,3,Eigen::Dynamic> Matrix3Xd;
	
    fsi_strip(lexer*,int);
	virtual ~fsi_strip();
	void start(lexer*,fdm*,ghostcell*,double);
	void initialize(lexer*,fdm*,ghostcell*,turbulence*);
    
    // Velocity interpolation: local contribution only, MPI reduction is done
    // once for all strips in fsi_strips::forcing (pack/unpack below)
    void interpolate_vel(lexer*,fdm*,ghostcell*,field&,field&,field&);
    int  numLagrangePoints() const;
    void pack_vel(double*) const;
    void unpack_vel(const double*);

    void update_points();
    void coupling_vel();
    void coupling_force(lexer*,double);
    // Spreads forces onto fx,fy,fz and accumulates the RANS forcing into eps0.
    // Ghost-cell exchange and epsget are done once for all strips by the caller.
    void distribute_forces(lexer*,fdm*,ghostcell*,field&,field&,field&,field&);
    void store_variables(lexer*);
    void print_ini(lexer *p);
    void print_stl(lexer*,fdm*,ghostcell*);
    void print_parameter(lexer*,fdm*,ghostcell*);
    
    void setFieldBC(Matrix3Xd&, Matrix3Xd&, Matrix4Xd&, Matrix4Xd&, Matrix4Xd&, Matrix3Xd&, Matrix4Xd&, Matrix3Xd&, double, int);
    void setConstantLoads(Matrix3Xd&, Matrix4Xd&, const Matrix3Xd&, const Matrix3Xd&, const Matrix4Xd&, const Matrix4Xd&);
    void setVariableLoads(Matrix3Xd&, Matrix4Xd&, const Matrix3Xd&, const Matrix3Xd&, const Matrix4Xd&, const Matrix4Xd&, const double);


private:

    void ini_parallel(lexer*,fdm*,ghostcell*);
    void get_cellsize(lexer*,fdm*,ghostcell*);
    void build_strip();
    double kernel_roma(const double&);
    bool near_subdomain(const Eigen::Vector3d&) const;
    void precompute_element_constants();
    void precompute_fluid_momentum();

    // Parallelisation
    int nstrip;
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;

    // Strip
    int Ne;
    vector<Matrix3Xd> Xil, Xil_0, lagrangePoints, lagrangeVel, lagrangeVelCoup, lagrangeForceCoup; 
    Matrix3Xd x_el, xdot_el, F_el, P_el, P_el_n, M_el, I_el, I_el_n;
    Matrix4Xd q_el, qdot_el;
    vector<Eigen::VectorXd> lagrangeArea;
    double rho_s,rho_f,dx_body,t_strip,t_strip_n,l_el,A_el,W_el;
    Eigen::Vector3d gravity_vec;
    bool thinStrip;

    // Cached quantities for setVariableLoads (previously recomputed in every
    // RHS evaluation of the beam solver, i.e. ~n_dim times per FD Jacobian)
    Eigen::VectorXd m_el_c;              // element fluid mass      (constant)
    Matrix3Xd s0_c;                      // static moment           (constant)
    vector<Eigen::Matrix3d> J0_c;        // quaternionic inertia    (constant)
    Matrix3Xd P_star_c, I_star_c;        // fluid momentum          (constant during one Integrate call)

    // Subdomain bounds extended by the kernel footprint (force spreading)
    double xlo_ext, xhi_ext, ylo_ext, yhi_ext, zlo_ext, zhi_ext;

    // Print
    int printcount_fsi;
	double printtime;
    Matrix3Xd tri_x, tri_y, tri_z;
    
    double starttime, endtime;
    
    turbulence *pturb;
};

#endif
