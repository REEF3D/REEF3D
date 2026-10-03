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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#ifndef SIXDOF_OBJ_H_
#define SIXDOF_OBJ_H_

#include"ddweno_f_nug.h"
#include<functional>
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field4a.h"
#include"field5.h"
#include"fieldint5.h"
#include"slice4.h"
#include"sliceint5.h"
#include"vtp3D.h"
#include"geo_raycast.h"
#include"6DOF_pto.h"
#include"6DOF_pto_joint.h"
#include"6DOF_rigidbody.h"
#include"6DOF_geometry.h"
#include<fstream>
#include<iostream>
#include<vector>
#include <Eigen/Dense>

class lexer;
class fdm;
class fdm2D;
class fdm_nhf;
class fdm_fnpf;
class ghostcell;
class reinidisc;
class nhflow_reinidisc_fsf;
class mooring;
class net_interface;
class sixdof_motionext;
 
using namespace std;

class sixdof_obj : public ddweno_f_nug, private vtp3D
{
public:
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	
    sixdof_obj(lexer*, ghostcell*, int);
	virtual ~sixdof_obj();
    
    // rigid-body core: state, kinematics and time integration (solver independent)
    sixdof_rigidbody rb;
    
    // surface geometry: hull triangles and pose transformation (solver independent)
    sixdof_geometry geom;
	
	void solve_eqmotion_cfd(lexer*,fdm*,ghostcell*,int,bool);
    
	void initialize_cfd(lexer*,fdm*,ghostcell*);
    void initialize_nhflow(lexer*,fdm_nhf*,ghostcell*);
    void initialize_shipwave(lexer*,ghostcell*,slice&,slice&);
    void initialize_wavemaker(lexer*,fdm_nhf*,ghostcell*,slice&,slice&);
    
	// Additional functions
    void transform(lexer*, fdm*, ghostcell*, bool);
    void update_forcing(lexer*, fdm*, ghostcell*,field&,field&,field&,field&,field&,field&,int);
    void hydrodynamic_forces_cfd(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void hydrodynamic_forces_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&,bool);
	
    void quat_matrices(lexer*);
    void update_position_3D(lexer*, fdm*, ghostcell*, bool);
    void update_position_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&, bool);
    void update_wavemaker_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&, bool);
    void update_position_2D(lexer*, ghostcell*,slice&);
    
    void solve_eqmotion_oneway_onestep(lexer*,ghostcell*,bool);
    
    // NHFLOW
    void solve_eqmotion_nhflow(lexer*,fdm_nhf*,ghostcell*,int,bool);
    void solve_eqmotion_oneway_nhflow(lexer*,ghostcell*,int,bool);
    void update_forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, slice&, int);
    void update_forcing_nhflow_wavemaker(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, slice&, int);
    void hydrodynamic_forces_nhflow_volume(lexer*, fdm_nhf*, ghostcell*,
                                           double*, double*, double*, slice&, int, bool);
    double Hsolidface_nhflow(lexer*, fdm_nhf*, int,int,int);
    
    // impermeable membranes (X 330) attached to the body: kinematics and forcing before the projection, loads after it
    void membrane_forcing_nhflow(lexer*,fdm_nhf*,ghostcell*,double,double*,double*,double*,slice&);
    void membrane_reaction_nhflow(lexer*,fdm_nhf*,ghostcell*,double,slice&,bool);
    bool membrane_iterated();
    void membrane_reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,double,double*,double*,double*,slice&);
    bool membrane_couple_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double,slice&,int);
    void membrane_stabilisation(lexer*,int);
    Eigen::Vector3d umem_n_=Eigen::Vector3d::Zero(), amem_n_=Eigen::Vector3d::Zero();
    double tmem_n_=-1.0;
    Eigen::Vector3d wmem_n_=Eigen::Vector3d::Zero(), almem_n_=Eigen::Vector3d::Zero();
    
    // porous floating body (X 16)
    void update_forcing_nhflow_porous(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, int);
    void porosity_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void porous_damping_nhflow(lexer*, int);
                         
                         
    double Sfx_n,Sfy_n,Sfz_n,SKx_n,SKy_n,SKz_n;
    int fictmass_ini;
    
    // print
    void saveTimeStep(lexer*,int);
    void print_parameter(lexer*,ghostcell*);
    void print_ini_vtp(lexer*,ghostcell*);
	void print_vtp(lexer*,ghostcell*);
    void print_normals_vtp(lexer*,ghostcell*);
    void print_ini_stl(lexer*,ghostcell*);
	void print_stl(lexer*,ghostcell*);
	void update_fbvel(lexer*,ghostcell*);
    
    // SFLOW
    double Hsolidface_2D(lexer*, int,int);
    void updateForcing_box(lexer*, ghostcell*, slice&);
    void updateForcing_stl(lexer*, ghostcell*, slice&, slice&);
    void updateForcing_oned(lexer*, ghostcell*, slice&);
    
    void update_forcing_sflow(lexer*, ghostcell*, slice&, slice&, slice&, slice&, slice&, slice&, int);
    
    void solve_eqmotion_sflow(lexer*,ghostcell*,int,bool);
    void solve_eqmotion_oneway_sflow(lexer*,ghostcell*,int,bool);
    
    double &Mass_fb;   // = rb.mass
    double Vfb, Rfb;
    
    // FNPF: resolved bodies in the sigma grid (6DOF_obj_fnpf*.cpp)
    void initialize_fnpf(lexer*, fdm_fnpf*, ghostcell*);
    void solve_eqmotion_fnpf(lexer*, ghostcell*, int, bool);
    void update_position_fnpf(lexer*, ghostcell*, bool);
    void print_fnpf(lexer*, ghostcell*, int);
    void ray_cast_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, slice&);
    void face_data_fnpf(lexer*, fdm_fnpf*, ghostcell*, int, double*, double*, double*, double*);
    void forces_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, double**, bool);
    // forces_fnpf in two parts for several grids (FNPF mesh refinement): the hull triangles
    // whose centroid own(x,y) accepts are integrated on grid (p,c) with the sampling distance
    // del, then the sums are reduced and stored
    struct fnpf_force_sum { double F[3], Mo[3], Am[36], Atot; };
    void forces_fnpf_zero(lexer*, fnpf_force_sum&);
    void forces_fnpf_sum(lexer*, fdm_fnpf*, double*, double**, bool, double, const std::function<bool(double,double)>*, fnpf_force_sum&);
    void forces_fnpf_set(lexer*, ghostcell*, fnpf_force_sum&, bool);
    double fnpf_dsm() const {return DSM;}
    bool fnpf_fixed(lexer*);
    void print_force_fnpf(lexer*);

    // read-only access for the SFLOW mesh refinement (sflow_amr)
    int amr_tricount() const {return tricount;}
    double **amr_tri(int d) {return d==0?tri_x:(d==1?tri_y:tri_z);}
    double amr_u(int n) const {return u_fb(n);}
    double amr_c(int n) const {return c_(n);}
    double amr_ramp_draft(lexer *p) {return ramp_draft(p);}
    slice& amr_fs() {return fs;}

    // NHFLOW mesh refinement (nhflow_amr): the body on a refined grid (level set FB of the grid
    // with its own ray-cast workspace IO, CL, CR of the grid size and the band spacing dsm)
    void ray_cast_nhflow_grid(lexer*, fdm_nhf*, ghostcell*, int*, int*, int*, double);
    double nhflow_dsm() const {return DSM;}
    // the grid that samples the load of a hull triangle with centroid (x,y): lexer, fdm and
    // water level (empty: the grid of the call)
    struct nhflow_grid { lexer *p; fdm_nhf *d; slice *WL; };
    std::function<nhflow_grid(double,double)> amr_grid_nhflow;
    // the hull triangles (X 185) are sized for the horizontal spacing times amr_hfac: the finest
    // level of a refinement zone around the body (set before the initialisation)
    double &amr_hfac;   // = geom.amr_hfac

private:

	void ini_parameter_stl(lexer*, fdm*, ghostcell*);
    void ini_fbvel(lexer*, ghostcell*);
    void maxvel(lexer*, ghostcell*);
    
    void externalForces_cfd(lexer*, fdm*, ghostcell*, double, bool);
    void externalForces_nhflow(lexer*, fdm_nhf*, ghostcell*, double, bool);
    void mooringForces(lexer*,  ghostcell*, double);
    void netForces_cfd(lexer*, fdm*, ghostcell*, double, bool);
    void netForces_nhflow(lexer*, fdm_nhf*, ghostcell*, double, bool);
    void update_forces(lexer*);
    
    double ramp_vel(lexer*);
    double ramp_draft(lexer*);
    
    void objects_create(lexer*, ghostcell*);
    void piston(lexer*, ghostcell*,int);
    void flap(lexer*, ghostcell*,int);
    void flap_double(lexer*, ghostcell*,int);
   
    void ini_parallel(lexer*, ghostcell*);
    
    double Hsolidface(lexer*, fdm*, int,int,int);
	double Hsolidface_t(lexer*, fdm*, int,int,int);
    
	
	void geometry_parameters(lexer*, fdm*, ghostcell*);
    void geometry_parameters_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void geometry_parameters_2D(lexer*, ghostcell*);
    void geometry_stl(lexer*, ghostcell*);
	void geometry_f(double&,double&,double&,double&,double&,double&,double&,double&,double&);
    void geometry_ls(lexer*, fdm*, ghostcell*);
    void geometry_ls_nhflow(lexer*, fdm_nhf*, ghostcell*);
    
    // force
    void forces_stl(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void forces_lsm(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void triangulation(lexer*, fdm*, ghostcell*, field&);
	void reconstruct(lexer*, fdm*, field&);
    void addpoint(lexer*,fdm*,int,int);
    void forces_lsm_calc(lexer* p, fdm *a, ghostcell *pgc,int,bool);
    
    void print_force(lexer*,fdm*,ghostcell*);
    void print_ini(lexer*,fdm*,ghostcell*);
    void print_vtp(lexer*,fdm*,ghostcell*);
    void pvtp(lexer*,int);
    
    void iniPosition_RBM(lexer*, ghostcell*);
    void update_trimesh_3D(lexer*, fdm*, ghostcell*, bool);
    void update_trimesh_nhflow(lexer*, fdm_nhf*, ghostcell*, bool);
    void update_trimesh_2D(lexer*, ghostcell*);
    void motionext_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&);
    void motionext_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&);
    

    // right-hand side of the rigid-body equations incl. prescribed motions (sixdof_motionext)
    void get_trans(lexer*, ghostcell*);
    void get_rot(lexer*);
    
    
    void rk2(lexer*, ghostcell*,int);
    void rk3(lexer*, ghostcell*,int);
    void rkls3(lexer*, ghostcell*,int);

   
   // ray cast 3D: geometry core kernels
    void ray_cast(lexer*, fdm*, ghostcell*);
    geo_raycast georay;
    void reini_RK2(lexer*, fdm*, ghostcell*, field&);
    
    // Raycast 3D
    fieldint5 cutl,cutr,fbio;
    // surface triangles: references into geom
    double **&tri_x,**&tri_y,**&tri_z,**&tri_x0,**&tri_y0,**&tri_z0;
    int &entity_sum;
    int *&tstart,*&tend;
    double xs,xe,ys,ye,zs,ze;
    int count, rayiter;
    double epsifb;
    const double epsi; 
    
    reinidisc *prdisc;
    net_interface *pnetinter;
    
	field4a f, frk1, L, dt; 
    int reiniter;
    
    
    // ray cast NHFLOW: geometry core kernels
    void ray_cast(lexer*, fdm_nhf*, ghostcell*);
    int  clip_facet_poly(lexer*,double,double,double,double,double,double,double,double,double,
                         double,double*,double*,double*);
    
    double zmin,zmax;
    double NB;
    
    // Reini NHFLOW
    void nhflow_reini_RK2(lexer*, fdm_nhf*, ghostcell*, double*);
    nhflow_reinidisc_fsf *pnhfrdisc;
    int *IO,*CR,*CL;
    double *FRK1,*DTT,*LL;
    
    // Forces NHFLOW
    void allocate(lexer*,fdm_nhf*,ghostcell*);
    void deallocate(lexer*,fdm_nhf*,ghostcell*);
    void print_force(lexer*,fdm_nhf*,ghostcell*);
    int *vert, *nflag;
    double *fsf;
    double A,Ax,Ay,Az;
    double Fx,Fy,Fz;
    
    double xp1,xp2,yp1,yp2,zp1,zp2;
    double sgnx,sgny,sgnz;
    double xloc,yloc,zloc;
	double xlocvel,ylocvel,zlocvel;
    double etaval,hspval;
    
    // -----
    // ray cast 2D
    void ray_cast_2D(lexer*, ghostcell*);
	void ray_cast_2D_io_x(lexer*, ghostcell*,int,int);
	void ray_cast_2D_io_ycorr(lexer*, ghostcell*,int,int);
    void ray_cast_2D_x(lexer*, ghostcell*,int,int);
	void ray_cast_2D_y(lexer*, ghostcell*,int,int);
    void ray_cast_2D_z(lexer*, ghostcell*,int,int);
    void reini_2D(lexer*,ghostcell*,slice&);
    void disc_2D(lexer*,ghostcell*,slice&);
    void time_preproc_2D(lexer*);
    
    slice4 press,lrk1,lrk2,K,dts,fs,Ls,Bs,Rxmin,Rxmax,Rymin,Rymax,draft;
    sliceint5 cl,cr,fsio;

    // Force NHFLOW
    void forces_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void force_calc_stl(lexer*, fdm_nhf*, ghostcell*, slice&,bool);
    void force_calc_stl2(lexer*, fdm_nhf*, ghostcell*, slice&,bool);
    void hydrodynamic_viscous_forces_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&,
            double&,double&,double&,double,double,double,double,double,double,double);
    void force_calc_lsm(lexer*, fdm_nhf*, ghostcell*,slice&);
    void triangulation(lexer*, fdm_nhf*, ghostcell*);
	void reconstruct(lexer*, fdm_nhf*);
	void addpoint(lexer*,fdm_nhf*,int,int);
	void finalize(lexer*,fdm_nhf*);
    double triangle_area(lexer*,double,double,double,double,double,double,double,double,double);
    double clip_edge(double,double);
    double clip_edge_vol(double,double,double);
    bool clip_facet(lexer*,double,double,double,double,double,double,double,double,double,
                double,double,double,double&,double&,double&,double&);
    void buoyancy_nhflow(lexer*, fdm_nhf*, ghostcell*, double,
                         double&, double&, double&, double&);
    
    // -----
    
    /* Rigid body state: references into the core (rb), so that the geometry, forcing and force
       routines keep their names
        - e: quaternions
        - h: angular momentum in body-fixed coordinates
        - c: position of mass centre in inertial system
        - p: linear momentum in inertial system
    */
    Eigen::Vector3d &p_, &c_, &h_, &dc_;
    Eigen::Vector4d &e_;
    Eigen::Matrix3d &R_, &I_, &quatRotMat;
    Eigen::Vector3d &omega_B, &omega_I;
    double &phi, &theta, &psi;
    
    Eigen::Matrix<double, 6, 1> u_fb;
    
    int &tricount, &entity_count;   // = geom
    
    double Uext, Vext, Wext, Pext, Qext, Rext;
    
    
    // extmotion
    sixdof_motionext *pmotion;
    
    // forces
    int **tri, **facet, *confac, *numfac,*numpt;
	double **ccpt, **pt, *ls;
	double   dV1,dV2,C1,C2,mi;
	int numtri,numvert, numtri_mem, numvert_mem;
	int countM,n,nn;
	int ccptcount,facount,check;
	int polygon_sum,polygon_num,vertice_num;
	const double zero,interfac;
    double eps;
    
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double xc,yc,zc;
    double nx,ny,nz,norm;
    double nxs,nys,nzs;
    double uval,vval,wval,pval,viscosity,density,phival;
    double du,dv,dw;
    double at,bt,ct,st;
    char name[100];
	
	fieldint5 vertice, nodeflag;
    field5 eta;
    
    int triangle_token,printnormal_count;
    
    double alpha[3],gamma[3],zeta[3];
    
    
    // Parallel	
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
   
    double kernel(const double&);

    // Print
    double curr_time;
    double printtime,printtimenormal;
    double *printtime_wT;
    int nCorr;
    int q,iin;
    float ffn;
    int offset[100];
    ofstream printpos,printforce,printvel;
    
    // Forces
    double Xext, Yext, Zext, Kext, Mext, Next;
    Eigen::Vector3d &Ffb_, &Mfb_;   // = rb.F, rb.M
    double Xe, Ye, Ze, Ke, Me, Ne;
    
    // porous floating body: drag reaction of the fluid on the skeleton (X 16)
    double Xd, Yd, Zd, Kd, Md, Nd;
    double Apor_fb, Bpor_fb;
    double Dpor_t, Dpor_r[3];

    // Mooring
	vector<double> X311_xen, X311_yen, X311_zen;
	vector<mooring*> pmooring;
	vector<double> Xme, Yme, Zme, Kme, Mme, Nme;    
    
    // Net Forces
    vector<double> Xne, Yne, Zne, Kne, Mne, Nne;   

    // Number
    int n6DOF;
    
    // Wavemaker
    double xwm1,zwm1,xwm2,zwm2;
    double *uwm,*wwm;
    
    void read_format_piston(lexer*,ghostcell*);
    void read_format_flap(lexer*,ghostcell*);
    void read_format_flap_double(lexer*,ghostcell*);
    
    double ts,te;
    double f0;
    int timecount,timecount_old;
    int rowcount,colcount;
    int colnum;
    int ptnum;
    double **kinematics;
    
    double DSM;
    
    // FNPF: classical RK4 stage derivatives and the added-mass coupling
    void rk4(lexer*, ghostcell*, int);
    void externalForces_fnpf(lexer*, ghostcell*, int, bool);
    void apply_added_mass(lexer*);
    bool p_fixed_dof(lexer*, int);
    Eigen::Matrix<double, 6, 6> Aadd_;
    bool am_on_ = false;

    // FNPF: power take-off (X 500 - X 504, 6DOF_obj_pto.cpp)
    void ini_pto(lexer*, ghostcell*);
    void pto_forces(lexer*, ghostcell*, int);
    void pto_implicit(lexer*, Eigen::Matrix<double,6,6>&, Eigen::Matrix<double,6,1>&);
    void pto_stability_check(lexer*);
    pto_joint_prismatic pto_joint_;
    pto_composite pto_;
    pto_state pto_s_;
    pto_output pto_out_;
    bool pto_on_ = false, pto_implicit_ = false, pto_warned_ = false;
    double pto_E_ = 0.0, pto_Pn_ = 0.0, pto_tn_ = -1.0;
    ofstream printpto;
};

#endif
