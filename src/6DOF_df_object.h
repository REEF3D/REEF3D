/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"gradient.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field5.h"
#include"fieldint5.h"
#include"vec.h"
#include<fstream>
#include<iostream>
#include<vector>
#include <Eigen/Dense>

class lexer;
class fdm;
class ghostcell;
class reinidisc;
class mooring;
class net;
class vrans;
 
using namespace std;

#ifndef SIXDOF_DF_OBJECT_H_
#define SIXDOF_DF_OBJECT_H_

class sixdof_df_object : public gradient
{
public:
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	
    sixdof_df_object(lexer*, fdm*, ghostcell*, int);
	virtual ~sixdof_df_object();
	
	virtual void solve_eqmotion(lexer*,fdm*,ghostcell*,int,vrans*,vector<net*>&);
	virtual void initialize(lexer*,fdm*,ghostcell*,vector<net*>&);
    
	// Additional functions
    void transform(lexer*, fdm*, ghostcell*, bool);
    void updateForcing(lexer*, fdm*, ghostcell*,double,field&,field&,field&,field&,field&,field&);
	void forces_stl(lexer*, fdm*, ghostcell*,double,field&,field&,field&);
    
    void saveTimeStep(lexer*,double);
    void print_parameter(lexer*,fdm*,ghostcell*);
    void print_ini_vtp(lexer*,fdm*,ghostcell*);
	void print_vtp(lexer*,fdm*,ghostcell*);
    void print_normals_vtp(lexer*,fdm*,ghostcell*);
    void print_ini_stl(lexer*,fdm*,ghostcell*);
	void print_stl(lexer*,fdm*,ghostcell*);
	void interface(lexer*, bool);
    
    double Mass_fb, Vfb, Rfb;

private:

	void ini_parameter_stl(lexer*, fdm*, ghostcell*);
    void ini_fbvel(lexer*, fdm*, ghostcell*);
    void maxvel(lexer*, fdm*, ghostcell*);
    
    void externalForces(lexer*, fdm*, ghostcell*, double, vrans*, vector<net*>&);
    void mooringForces(lexer*, fdm*, ghostcell*, double);
    void netForces(lexer*, fdm*, ghostcell*, double, vrans*, vector<net*>&);
    void updateForces(fdm*);
    
    void objects_create(lexer*, fdm*, ghostcell*);
    void objects_allocate(lexer*, fdm*, ghostcell*);
	void geometry_refinement(lexer*,ghostcell*);
	void create_triangle(double&,double&,double&,double&,double&,double&,double&,double&,double&,const double&,const double&,const double&);
	void box(lexer*, fdm*, ghostcell*,int);
	void cylinder_x(lexer*, fdm*, ghostcell*,int);
	void cylinder_y(lexer*, fdm*, ghostcell*,int);
	void cylinder_z(lexer*, fdm*, ghostcell*,int);
	void wedge_sym(lexer*, fdm*, ghostcell*,int);
    void wedge(lexer*, fdm*, ghostcell*,int);
    void hexahedron(lexer*, fdm*, ghostcell*,int);
    void read_stl(lexer*, fdm*, ghostcell*);
    void triangle_switch_lsm(lexer*, fdm*, ghostcell*);
    void triangle_switch_ray(lexer*, fdm*, ghostcell*);
   
    void ini_parallel(lexer*, fdm*, ghostcell*);
    
    double Hsolidface(lexer*, fdm*, int,int,int);
	double Hsolidface_t(lexer*, fdm*, int,int,int);
	
	void geometry(lexer*, fdm*, ghostcell*);
    void geometry_stl(lexer*, fdm*, ghostcell*);
	void geometry_f(double&,double&,double&,double&,double&,double&,double&,double&,double&);
    void geometry_ls(lexer*, fdm*, ghostcell*);
    
    void iniPosition_RBM(lexer*, fdm*, ghostcell*);
    void update_Position(lexer*, fdm*, ghostcell*, bool);
    void prescribedMotion(lexer*, fdm*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&);
    void quat_matrices(const Eigen::Vector4d&);

    void get_trans(lexer*, fdm*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
    void get_rot(Eigen::Vector3d&, Eigen::Vector4d&, const Eigen::Vector3d&, const Eigen::Vector4d&);
    Eigen::Matrix3d quatRotMat;
    
    /*Eigen::Vector3d ab4_3
    (   
        lexer*,
        const Eigen::Vector3d&, 
        const Eigen::Vector3d&, 
        const Eigen::Vector3d&, 
        const Eigen::Vector3d&,
        double alpha
    );
    Eigen::Vector4d ab4_4
    (
        lexer*,
        const Eigen::Vector4d&, 
        const Eigen::Vector4d&, 
        const Eigen::Vector4d&, 
        const Eigen::Vector4d&,
        double alpha
    );
    Eigen::Vector3d am4_3
    (
        lexer*,
        const Eigen::Vector3d&, 
        const Eigen::Vector3d&, 
        const Eigen::Vector3d&, 
        const Eigen::Vector3d&,
        double alpha
    );
    Eigen::Vector4d am4_4
    (
        lexer*,
        const Eigen::Vector4d&, 
        const Eigen::Vector4d&, 
        const Eigen::Vector4d&, 
        const Eigen::Vector4d&,
        double alpha
    );
    void abam4(lexer*, fdm*, ghostcell*, double);
    void rk2(lexer*, fdm*, ghostcell*,double);
    void rk4(lexer*, fdm*, ghostcell*,double);
    */
    
    void rk2(lexer*, fdm*, ghostcell*,int);
    void rk3(lexer*, fdm*, ghostcell*,int);
    void rkls3(lexer*, fdm*, ghostcell*,int);

    void rotation_tri(lexer*,double,double,double,double&,double&,double&, const double&, const double&, const double&);
   
    void ray_cast(lexer*, fdm*, ghostcell*);
	void ray_cast_io_x(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_io_ycorr(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_io_zcorr(lexer*, fdm*, ghostcell*,int,int);
    void ray_cast_x(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_y(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_z(lexer*, fdm*, ghostcell*,int,int);
    void reini_AB2(lexer*, fdm*, ghostcell*, field&);
    
    

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
    
    int tricount, entity_count;
    
    double phi, theta, psi;
    
    double Uext, Vext, Wext, Pext, Qext, Rext;
    
    double dtn1, dtn2, dtn3;

    
    // Raycast
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
    
	// Parallel	
	double *xstart, *xend, *ystart, *yend, *zstart, *zend;
	
    // Reini
    reinidisc *prdisc;
	vec f, frk1, L, dt; 
    int reiniter;
   
    double kernel(const double&);

    // Print
    double curr_time;
    double printtime,printtimenormal;
    double *printtime_wT;
    int nCorr;
    int q,iin;
    float ffn;
    int offset[100];
    
    
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
