/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#ifndef SIXDOF_OBJ_H_
#define SIXDOF_OBJ_H_

#include"ddweno_f_nug.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field4a.h"
#include"field5.h"
#include"fieldint5.h"
#include"slice4.h"
#include"sliceint5.h"
#include<fstream>
#include<iostream>
#include<vector>
#include <Eigen/Dense>

class lexer;
class fdm;
class fdm2D;
class fdm_nhf;
class ghostcell;
class reinidisc;
class nhflow_reinidisc_fsf;
class mooring;
class net;
class vrans;
class sixdof_motionext;
 
using namespace std;

class sixdof_obj : public ddweno_f_nug
{
public:
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	
    sixdof_obj(lexer*, ghostcell*, int);
	virtual ~sixdof_obj();
	
	virtual void solve_eqmotion(lexer*,fdm*,ghostcell*,int,vrans*,vector<net*>&);
    
	void initialize_cfd(lexer*,fdm*,ghostcell*,vector<net*>&);
    void initialize_nhflow(lexer*,fdm_nhf*,ghostcell*,vector<net*>&);
    void initialize_shipwave(lexer*,ghostcell*,slice&,slice&);
    
	// Additional functions
    void transform(lexer*, fdm*, ghostcell*, bool);
    void update_forcing(lexer*, fdm*, ghostcell*,field&,field&,field&,field&,field&,field&,int);
    void hydrodynamic_forces_cfd(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void hydrodynamic_forces_nhflow(lexer*, fdm_nhf*, ghostcell*,bool);
	
    void quat_matrices(lexer*);
    void update_position_3D(lexer*, fdm*, ghostcell*, bool);
    void update_position_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&, bool);
    void update_position_2D(lexer*, ghostcell*,slice&);
    
    void solve_eqmotion_oneway_onestep(lexer*,ghostcell*);
    
    // NHFLOW
    virtual void solve_eqmotion_nhflow(lexer*,fdm_nhf*,ghostcell*,int,vrans*,vector<net*>&);
    void solve_eqmotion_oneway_nhflow(lexer*,ghostcell*,int);
    void update_forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, slice&, int);
    
    double Hsolidface_nhflow(lexer*, fdm_nhf*, int,int,int);
    
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
    
    void solve_eqmotion_sflow(lexer*,ghostcell*,int);
    void solve_eqmotion_oneway_sflow(lexer*,ghostcell*,int);
    
    double Mass_fb, Vfb, Rfb;

private:

	void ini_parameter_stl(lexer*, fdm*, ghostcell*);
    void ini_fbvel(lexer*, ghostcell*);
    void maxvel(lexer*, ghostcell*);
    
    void externalForces(lexer*, fdm*, ghostcell*, double, vrans*, vector<net*>&);
    void externalForces_nhflow(lexer*, fdm_nhf*, ghostcell*, double, vrans*, vector<net*>&);
    void mooringForces(lexer*,  ghostcell*, double);
    void netForces(lexer*, fdm*, ghostcell*, double, vrans*, vector<net*>&);
    void update_forces(lexer*);
    
    double ramp_vel(lexer*);
    double ramp_draft(lexer*);
    
    void objects_create(lexer*, ghostcell*);
    void objects_allocate(lexer*, ghostcell*);
	void geometry_refinement(lexer*,ghostcell*);
	void create_triangle(double&,double&,double&,double&,double&,double&,double&,double&,double&,const double&,const double&,const double&);
	void box(lexer*, ghostcell*,int);
	void cylinder_x(lexer*, ghostcell*,int);
	void cylinder_y(lexer*, ghostcell*,int);
	void cylinder_z(lexer*, ghostcell*,int);
	void wedge_sym(lexer*, ghostcell*,int);
    void wedge(lexer*, ghostcell*,int);
    void hexahedron(lexer*, ghostcell*,int);
    void read_stl(lexer*, ghostcell*);
    void triangle_switch_lsm(lexer*, ghostcell*);
    void triangle_switch_ray(lexer*, ghostcell*);
   
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
    void pvtp(lexer*,fdm*,ghostcell*);
    void header(lexer*,fdm*,ghostcell*);
    void name_iter(lexer*,fdm*,ghostcell*);
    void name_time(lexer*,fdm*,ghostcell*);
    void piecename(lexer*,fdm*,ghostcell*,int);
    
    
    void iniPosition_RBM(lexer*, ghostcell*);
    void update_Euler_angles(lexer*, ghostcell*);
    void update_trimesh_3D(lexer*, fdm*, ghostcell*, bool);
    void update_trimesh_nhflow(lexer*, fdm_nhf*, ghostcell*, bool);
    void update_trimesh_2D(lexer*, ghostcell*);
    void motionext_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&);
    void motionext_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&);
    

    void get_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&);
    void get_rot(lexer*,Eigen::Vector3d&, Eigen::Vector4d&, Eigen::Vector3d&, Eigen::Vector4d&);
    Eigen::Matrix3d quatRotMat;
    
    
    void rk2(lexer*, ghostcell*,int);
    void rk3(lexer*, ghostcell*,int);
    void rkls3(lexer*, ghostcell*,int);

    void rotation_tri(lexer*,double,double,double,double&,double&,double&, const double&, const double&, const double&);
   
   // ray cast 3D
    void ray_cast(lexer*, fdm*, ghostcell*);
	void ray_cast_io_x(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_io_ycorr(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_io_zcorr(lexer*, fdm*, ghostcell*,int,int);
    void ray_cast_x(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_y(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_z(lexer*, fdm*, ghostcell*,int,int);
    void ray_cast_direct(lexer*, fdm*, ghostcell*,int,int);
    void reini_AB2(lexer*, fdm*, ghostcell*, field&);
    void reini_RK2(lexer*, fdm*, ghostcell*, field&);
    
    // Raycast 3D
    fieldint5 cutl,cutr,fbio;
    double **tri_x,**tri_y,**tri_z,**tri_x0,**tri_y0,**tri_z0;
    int *tri_switch,*tri_switch_id,*tri_switch_local,*tri_switch_local_id;
    int tricount_local,*tricount_local_list,*tricount_local_displ;
    int tricount_switch_total;
	vector<vector<double> > tri_x_r;
	vector<vector<double> > tri_y_r;
	vector<vector<double> > tri_z_r;
    double xs,xe,ys,ye,zs,ze;
    int entity_sum, count, rayiter;
    int *tstart,*tend;
    double epsifb;
    const double epsi; 
    
    reinidisc *prdisc;
	field4a f, frk1, L, dt; 
    int reiniter;
    
    // ray cast NHFLOW
    void ray_cast(lexer*, fdm_nhf*, ghostcell*);
    void ray_cast_io(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_x(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_y(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_z(lexer*, fdm_nhf*, ghostcell*,int,int);
    
    double zmin,zmax;
    
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
    void force_calc_stl(lexer*, fdm_nhf*, ghostcell*,bool);
    void force_calc_lsm(lexer*, fdm_nhf*, ghostcell*);
    void triangulation(lexer*, fdm_nhf*, ghostcell*);
	void reconstruct(lexer*, fdm_nhf*);
	void addpoint(lexer*,fdm_nhf*,int,int);
	void finalize(lexer*,fdm_nhf*);
    double triangle_area(lexer*,double,double,double,double,double,double,double,double,double);
    
    // -----
    
    /* Rigid body motion
        - e: quaternions
        - h: angular momentum in body-fixed coordinates
        - c: position of mass centre in inertial system
        - p: velocity of mass centre in inertial system
    */
    Eigen::Vector3d p_, pk_, pn1_, pn2_, pn3_, dp_, dpk_, dpn1_, dpn2_, dpn3_; 
    Eigen::Vector3d c_, ck_, cn1_, cn2_, cn3_, dc_, dck_, dcn1_, dcn2_, dcn3_;
    Eigen::Vector3d h_, hk_, hn1_, hn2_, hn3_, dh_, dhk_, dhn1_, dhn2_, dhn3_;
    Eigen::Vector4d e_, ek_, en1_, en2_, en3_, de_, dek_, den1_, den2_, den3_;
    Eigen::Matrix<double, 3, 4> E_, G_, Gdot_;
    
    Eigen::MatrixXd deltad_, delta_, deltan1_, deltan2_, deltan3_;

    Eigen::Matrix3d R_, I_, Rinv_;
    
    Eigen::Vector3d omega_B, omega_I;
    
    Eigen::Matrix<double, 6, 1> u_fb;
    
    int tricount, entity_count;
    
    double phi, theta, psi;
    
    double Uext, Vext, Wext, Pext, Qext, Rext;
    
    double dtn1, dtn2, dtn3;
    
    
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
    char name[100],pname[100];
	
	fieldint5 vertice, nodeflag;
    field5 eta;
    
    
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
    Eigen::Vector3d Ffb_, Mfb_;
    double Xe, Ye, Ze, Ke, Me, Ne;


    // Mooring
	vector<double> X311_xen, X311_yen, X311_zen;
	vector<mooring*> pmooring;
	vector<double> Xme, Yme, Zme, Kme, Mme, Nme;    
    
    // Net
    vector<double> Xne, Yne, Zne, Kne, Mne, Nne;   

    // Number
    int n6DOF;
    
    int triangle_token,printnormal_count;
    
    double alpha[3],gamma[3],zeta[3];
};

#endif
