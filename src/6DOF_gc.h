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
Authors: 
    Hans Bihs: Euler angle implementation
    Tobias Martin: Euler parameter implementation
--------------------------------------------------------------------*/

#include"6DOF.h"
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

class reinidisc;
class mooring;
class net;
 
using namespace std;

#ifndef SIXDOF_GC_H_
#define SIXDOF_GC_H_

class sixdof_gc : public sixdof, public gradient
{
public:
	sixdof_gc(lexer*, fdm*, ghostcell*);
	virtual ~sixdof_gc();
	
	virtual void start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&);
	virtual void initialize(lexer*,fdm*,ghostcell*,vector<net*>&);
    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);

private:
	void ini_parameter(lexer*, fdm*, ghostcell*);
	void interface(lexer*, bool);
	void start_Euler(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&);
	void start_Quaternion(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&);
	void reini_AB2(lexer*, fdm*, ghostcell*, field&);
	void position_ini(lexer*, fdm*, ghostcell*);
	void position_ini_quaternion(lexer*, fdm*, ghostcell*);
	void ray_cast(lexer*, fdm*, ghostcell*);
	void ray_cast_io_x(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_io_ycorr(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_io_zcorr(lexer*, fdm*, ghostcell*,int,int);
    void ray_cast_x(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_y(lexer*, fdm*, ghostcell*,int,int);
	void ray_cast_z(lexer*, fdm*, ghostcell*,int,int);
	void center_gravity(lexer*, fdm*, ghostcell*);
	void moments_inertia(lexer*, fdm*, ghostcell*);
	void geometry_ini(lexer*, fdm*, ghostcell*);
	void geometry_ini_f(double&,double&,double&,double&,double&,double&,double&,double&,double&);
	void fb_position(lexer*, fdm*, ghostcell*);
	void fb_position_quaternion(lexer*, fdm*, ghostcell*,const std::vector<double>&);
	void maxvel(lexer*, fdm*, ghostcell*);
	
	void objects(lexer*, fdm*, ghostcell*);
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
    void read_motionvec(lexer*, fdm*, ghostcell*);
	
    void rotation_stl(lexer*,double&,double&,double&);
    void rotation_stl_quaternion(lexer*,double,double,double,double&,double&,double&, const double&, const double&, const double&);
	
    void motion_ext(lexer*, fdm*, ghostcell*);
	void motion_ext_quaternion(lexer*, fdm*, ghostcell*);
	void motion_fixed(lexer*, fdm*, ghostcell*);
    void motion_vec(lexer*, fdm*, ghostcell*);
	void preventMotion(lexer*);
	
	void print_ini_vtp(lexer*,fdm*,ghostcell*);
	void print_vtp(lexer*,fdm*,ghostcell*);
    void print_ini_stl(lexer*,fdm*,ghostcell*);
	void print_stl(lexer*,fdm*,ghostcell*);
	void print_E_position(lexer*,fdm*,ghostcell*);
	void print_E_velocity(lexer*,fdm*,ghostcell*);
	void print_E_force(lexer*,fdm*,ghostcell*);
	void print_S_force(lexer*,fdm*,ghostcell*);
	void rotation(double&,double&,double&,double,double,double);
	void rotation_quaternion(double&,double&,double&,double,double,double);
	
	// forces
	void fluidForces(lexer*,fdm*,ghostcell*);
    void forces_lsm(lexer*,fdm*,ghostcell*);
	void forces_stl(lexer*,fdm*,ghostcell*);
    void forces_triang(lexer*,fdm*,ghostcell*);
    void forces_triang_triangulation(lexer*,fdm*,ghostcell*);
    void forces_triang_reconstruction(lexer*,fdm*,ghostcell*);
    void forces_triang_addpoint(lexer*,fdm*,int,int);
    void forces_triang_finalize(lexer*,fdm*,ghostcell*);
	void mooringForces(lexer*, fdm*, ghostcell*);
    void netForces(lexer*, fdm*, ghostcell*, double, vrans*, vector<net*>&);
	
	void solve(lexer*,fdm*,ghostcell*);
	void solve_quaternion();
    std::vector<double> get_R(const std::vector<double> &);
	std::vector<double> get_e(const std::vector<double> &, const std::vector<double> &);
	std::vector<double> get_h(const std::vector<double> &, const std::vector<double> &);
    std::vector<double> rotation_R(const std::vector<double>&);
	void update();
    void update_quaternion();
	void transform_vec_ES(double,double,double,double &,double &,double &);
	void transform_vec_SE(double,double,double,double &,double &,double &);
	void transform_angle_ES(double,double,double,double &,double &,double &);
	void transform_angle_SE(double,double,double,double &,double &,double &);
	
	
	void solidUpdate(lexer*,fdm*,ghostcell*,const std::vector<double>&);
	void forceUpdate(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&);
	
	// ------ 

	fieldint5 cutl,cutr,fbio;
    
	
	int conv(double);
	
	const double epsifb;
	
	double Vfb,Mfb,Rfb;
	double xg,yg,zg;
	double xg_s,yg_s,zg_s;
	double xg_sn,yg_sn,zg_sn;
	double xgn,ygn,zgn;
	double dxg,dyg,dzg;
	double Ix,Iy,Iz;
	double xorig,yorig,zorig;
    
    Eigen::Matrix3d quatRotMat;
	
    std::vector< std::vector<double> > R_;
    std::vector< std::vector<double> > I_;
    std::vector<double> e_;
    std::vector<double> en_;
    std::vector<double> enn_;
    std::vector<double> ennn_;
    std::vector<double> ennnn_;
    std::vector<double> trunc_;
    std::vector<double> L_;
    std::vector<double> Ln_;
    std::vector<double> Lnn_;
    
	double Ue,Ve,We;
	double Pe,Qe,Re;
	double dUe,dVe,dWe;
	double dPe,dQe,dRe;
	
	double Us,Vs,Ws;
	double Ps,Qs,Rs;
	double dUs,dVs,dWs;
	double dPs,dQs,dRs;
	
	double Uext,Vext,Wext;
	double Pext,Qext,Rext;
	
	field4 phin;
	
	double Usn,Vsn,Wsn;
	double Psn,Qsn,Rsn;
	double Usnn,Vsnn,Wsnn;
	double Psnn,Qsnn,Rsnn;
	double Usnnn,Vsnnn,Wsnnn;
	double Psnnn,Qsnnn,Rsnnn;
	
	double dUen,dVen,dWen;
	double dUenn,dVenn,dWenn;
	double dUennn,dVennn,dWennn;	
	double Uen,Ven,Wen;
	double Pen,Qen,Ren;
	double Uenn,Venn,Wenn;
	double Penn,Qenn,Renn;
	double Uennn,Vennn,Wennn;
	double Pennn,Qennn,Rennn;
	
	double dUsn,dVsn,dWsn;
	double dPsn,dQsn,dRsn;
	double dUsnn,dVsnn,dWsnn;
	double dPsnn,dQsnn,dRsnn;
	double dUsnnn,dVsnnn,dWsnnn;
	double dPsnnn,dQsnnn,dRsnnn;


	double dxg_sum,dyg_sum,dzg_sum;
	double dphi_sum,dtheta_sum,dpsi_sum;
	double phi,theta,psi;
	double phi_s,theta_s,psi_s;
	double phi_sn,theta_sn,psi_sn;
	double phi_en,theta_en,psi_en;
	double dphi,dtheta,dpsi;
	double phi1,theta1,psi1;
	double phi2,theta2,psi2;
	double Xe,Ye,Ze;
	double Ke,Me,Ne;
	double Xs,Ys,Zs;
	double Ks,Ms,Ns;
	
	double H;
	double rx,ry,rz;
	double sx,sy,sz;
	double printtime,starttime,endtime;
    
	int istart,iend,jstart,jend,kstart,kend;
	int q,count;
	double xs,xe,ys,ye,zs,ze;
	int ts,te;
	int *tstart,*tend;
	int tricount,trisum;
    int ccptcount,numtri;

	// mooring
	vector<double> X311_xen, X311_yen, X311_zen;
	vector<mooring*> pmooring;
	vector<double> Xme, Yme, Zme, Kme, Mme, Nme;

    // net
	vector<net*> pnet;
    vector<double> Xne, Yne, Zne, Kne, Mne, Nne;    
    
    // triangulation
	double **tri_x,**tri_y,**tri_z;
	double **tri_xn,**tri_yn,**tri_zn;
	vector<vector<double> > tri_x_r;
	vector<vector<double> > tri_y_r;
	vector<vector<double> > tri_z_r;
    int **tri, **facet, *confac, *numfac,*numpt;
    int **ijk;
	double **ccpt, **pt, *ls;
    double ***ccell,**lscc;
	int *ccnode,**ccid,*ccflag,**vertice_cc,***ccijk;
    
    double zero;
    const double epsi;
    int check,facount,countM;
    int numvert,numtri_mem,numvert_mem,polygon_num;
    int nn;
    int entity_sum;
    int entity_count;
    
    fieldint5 vertice, nodeflag;
	field5 eta;
    
    // motion
    double **motion;
	double vecx,vecy,vecz;
    double nvecx,nvecy,nvecz;
    double ts_motion,te_motion;
    int tcount_motion;
    
	// reini	
	reinidisc *prdisc;
	vec f,frk1,frk2,L;
	
	vec dt;
	int reiniter;
	double xmin,xplus,ymin,yplus,zmin,zplus;
	double dstx,dsty,dstz,lsSig,dnorm,op,lsv,sign;
	
	ofstream eposout, evelout, eforceout, sforceout;
    
    int rayiter;

    void print_forces_vtp(lexer*,fdm*,ghostcell*);
    void forces_pvtp(lexer*,fdm*,ghostcell*);
    void forces_header(lexer*,fdm*,ghostcell*);
    void name_iter(lexer*,fdm*,ghostcell*);
    void name_time(lexer*,fdm*,ghostcell*);
    void piecename(lexer*,fdm*,ghostcell*, int);
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int forceprintcount;
};

#endif
