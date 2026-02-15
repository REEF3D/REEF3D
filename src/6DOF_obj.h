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
#include"field_header.h"
#include"fieldint4.h"
#include"slice4.h"
#include"sliceint5.h"
#include"vtp3D.h"
#include"geo_raycast.h"
#include"6DOF_pto.h"
#include"6DOF_pto_joint.h"
#include"6DOF_pto_controller.h"
#include"6DOF_rigidbody.h"
#include"6DOF_geometry.h"
#include"6DOF_load.h"
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

class sixdof_obj : public ddweno_f_nug, protected vtp3D
{
public:
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
    sixdof_obj(lexer*, ghostcell*, int);
	virtual ~sixdof_obj();
    
    // rigid-body core: state, kinematics and time integration (solver independent)
    sixdof_rigidbody rb;
    
    // mesh refinement with subcycling (G 7 1, 6DOF_obj_amr.cpp): the coarser levels step with a
    // predicted copy of the body; the state is saved before and put back after
    void amr_save();
    void amr_restore(lexer*, ghostcell*);
    sixdof_rigidbody rb_amr;
    
    // surface geometry: hull triangles and pose transformation (solver independent)
    sixdof_geometry geom;
    
    // external load models (e.g. the ship module, X 350), owned by the body, and the fluid access
    // of the solver coupling during update_forces (nullptr: none)
    vector<sixdof_load*> pload;
    sixdof_fluid *pfluid = nullptr;
	
    
    
	
    void quat_matrices(lexer*);
    
    void solve_eqmotion_oneway_onestep(lexer*,ghostcell*,bool);
    
    
    
                         
                         
    double Sfx_n,Sfy_n,Sfz_n,SKx_n,SKy_n,SKz_n;
    int fictmass_ini;
    
    // print
    void saveTimeStep(lexer*,int);
    void print_parameter(lexer*,ghostcell*);
    void print_ini_vtp(lexer*,ghostcell*);
	void print_vtp(lexer*,ghostcell*);
    void print_ini_stl(lexer*,ghostcell*);
	void print_stl(lexer*,ghostcell*);
	void update_fbvel(lexer*,ghostcell*);
    
    
    
    
    double &Mass_fb;   // = rb.mass
    double Vfb, Rfb;
    

    // read-only access for the SFLOW mesh refinement (sflow_amr)
    int amr_tricount() const {return tricount;}
    double **amr_tri(int d) {return d==0?tri_x:(d==1?tri_y:tri_z);}
    double amr_u(int n) const {return u_fb(n);}
    double amr_c(int n) const {return c_(n);}
    double amr_ramp_draft(lexer *p) {return ramp_draft(p);}

    // the hull triangles (X 185) are sized for the horizontal spacing times amr_hfac: the finest
    // level of a refinement zone around the body (set before the initialisation)
    double &amr_hfac;   // = geom.amr_hfac

protected:

    void ini_fbvel(lexer*, ghostcell*);
    void maxvel(lexer*, ghostcell*);
    
    void mooringForces(lexer*,  ghostcell*, double);
    void update_forces(lexer*);
    
    double ramp_vel(lexer*);
    double ramp_draft(lexer*);
    
    void objects_create(lexer*, ghostcell*);
    void piston(lexer*, ghostcell*,int);
    void flap(lexer*, ghostcell*,int);
    void flap_double(lexer*, ghostcell*,int);
   
    void ini_parallel(lexer*, ghostcell*);
    
    
	
    void geometry_stl(lexer*, ghostcell*);
    void geometry_parameters_surface(lexer*, ghostcell*, double);   // mass properties from the surface triangles
	void geometry_f(double&,double&,double&,double&,double&,double&,double&,double&,double&);
    
    // force
    
    
    void iniPosition_RBM(lexer*, ghostcell*);
    void motionext_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&);
    void motionext_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&);
    

    // right-hand side of the rigid-body equations incl. prescribed motions (sixdof_motionext)
    void get_trans(lexer*, ghostcell*);
    void get_rot(lexer*);
    
    
    void rk2(lexer*, ghostcell*,int);
    void rk3(lexer*, ghostcell*,int);
    void rkls3(lexer*, ghostcell*,int);

   
   // ray cast 3D: geometry core kernels
    geo_raycast georay;
    
    // Raycast 3D
    // surface triangles: references into geom
    double **&tri_x,**&tri_y,**&tri_z,**&tri_x0,**&tri_y0,**&tri_z0;
    int &entity_sum;
    int *&tstart,*&tend;
    double xs,xe,ys,ye,zs,ze;
    int count, rayiter;
    double epsifb;
    const double epsi; 
    
    net_interface *pnetinter;
    
    int reiniter;
    
    
    
    
    double zmin,zmax;
    double NB;
    
    
    double A,Ax,Ay,Az;
    double Fx,Fy,Fz;
    
    double xp1,xp2,yp1,yp2,zp1,zp2;
    double sgnx,sgny,sgnz;
    double xloc,yloc,zloc;
	double xlocvel,ylocvel,zlocvel;
    double etaval,hspval;
    
    // -----
    

    
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
    void apply_added_mass(lexer*, const Eigen::Matrix<double,6,6>* = nullptr);
    bool p_fixed_dof(lexer*, int);
    Eigen::Matrix<double, 6, 6> Aadd_;
    bool am_on_ = false;

    // FNPF: power take-off (X 500 - X 508, 6DOF_obj_pto.cpp)
    void ini_pto(lexer*, ghostcell*);
    void pto_forces(lexer*, ghostcell*, int);
    void pto_implicit(lexer*, Eigen::Matrix<double,6,6>&, Eigen::Matrix<double,6,1>&);
    void pto_stability_check(lexer*);
    pto_joint_prismatic pto_joint_;
    pto_composite pto_;
    std::vector<std::unique_ptr<pto_controller>> pto_ctrl_;
    void pto_tuned_passive(lexer*, ghostcell*);
    pto_state pto_s_;
    pto_output pto_out_;
    bool pto_on_ = false, pto_implicit_ = false, pto_warned_ = false;
    double pto_E_ = 0.0, pto_Estep_ = 0.0;
    ofstream printpto;
};

#endif
