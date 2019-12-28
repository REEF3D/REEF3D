/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include <mpi.h>
#include"norm_vec.h"
#include"boundarycheck.h"


class fdm;
class fdm2D;
class lexer;
class field;
class fieldint;
class slice;
class sliceint;
class vec;
class vec2D;
class cpt;
class cpt2D;
class density_f;
class reini;
class convection;
class ioflow;

#ifndef GHOSTCELL_H_
#define GHOSTCELL_H_

using namespace std;

class ghostcell : public norm_vec, public boundarycheck
{
public:
	ghostcell(int&,char**,lexer*);
	virtual ~ghostcell();

	void gcini(lexer* p);

	void start1(lexer*,field&, int);
	void start2(lexer*,field&, int);
	void start3(lexer*,field&, int);
	void start4(lexer*,field&, int);
	void start4a(lexer*,field&, int);

	void start4V(lexer*,vec&,int);
	void start4aV(lexer*,vec&,int);
    void start6V(lexer*,vec&,int);
    void start7V(lexer*,double*,sliceint&, int);


	void final();
	double globalsum(double);
	int globalisum(int);
	double globalmax(double);
	double globalmin(double);
	int globalimax(int);
	int globalimin(int);
	double timesync(double);
	void globalctrl(lexer*);
	void dgcpol(lexer*,field&, int**,int, int);
    
    void dgcpol1(lexer*,field&, int);
    void dgcpol2(lexer*,field&, int);
    void dgcpol3(lexer*,field&, int);
    void dgcpol4(lexer*,field&, int);
    void dgcpol4a(lexer*,field&, int);
    void dgcpol6(lexer*,field&, int);
    
    void dgcslpol(lexer*, slice&, int**,int, int);
    void dgcslpol1(lexer*, slice&);
    void dgcslpol2(lexer*, slice&);
    void dgcslpol4(lexer*, slice&);
	void parapls(lexer*,double**,double**,int*,int*);

//  Update
	void facenbx(lexer*, fieldint&, int*);
	void flagx(lexer*,int*);
    void flagx7(lexer*,int*);
    void sigmax(lexer*,fdm*,double*);
    void rangex(lexer*,int*,int);
	void gcxupdate(lexer*);

    void dgcslini1(lexer*);
	void dgcslini2(lexer*);
	void dgcslini3(lexer*);
	void dgcslini4(lexer*);

	void cval_update1(lexer*,fdm*,fieldint&);
	void cval_update2(lexer*,fdm*,fieldint&);
	void cval_update3(lexer*,fdm*,fieldint&);
	void cval_update4(lexer*,fdm*,fieldint&);
    void cval_update4a(lexer*,fdm*,fieldint&);
    void cval_update6(lexer*,fdm*,fieldint&);

    void rownum4_update(lexer*,fieldint&);
    void rownum7_update(lexer*,int*);


	void cval_gcb1(lexer*,fdm*,fieldint&);
    void cval_gcb2(lexer*,fdm*,fieldint&);
    void cval_gcb3(lexer*,fdm*,fieldint&);
    void cval_gcb4(lexer*,fdm*,fieldint&);
    void cval_gcb4a(lexer*,fdm*,fieldint&);
    void cval_gcb6(lexer*,fdm*,fieldint&);

	void cval_gcpara1(lexer*,fdm*,fieldint&);
    void cval_gcpara2(lexer*,fdm*,fieldint&);
    void cval_gcpara3(lexer*,fdm*,fieldint&);
    void cval_gcpara4(lexer*,fdm*,fieldint&);
    void cval_gcpara4a(lexer*,fdm*,fieldint&);
    void cval_gcpara6(lexer*,fdm*,fieldint&);

    void column_pt1_update(lexer*,fdm*);
    void column_pt2_update(lexer*,fdm*);
    void column_pt3_update(lexer*,fdm*);
    void column_pt4_update(lexer*,fdm*);
    void column_pt4a_update(lexer*,fdm*);
    void column_pt6_update(lexer*,fdm*);

	void column_pt1(lexer*,fdm*,fieldint&);
	void column_pt2(lexer*,fdm*,fieldint&);
	void column_pt3(lexer*,fdm*,fieldint&);
	void column_pt4(lexer*,fdm*,fieldint&);
    void column_pt4a(lexer*,fdm*,fieldint&);
    void column_pt6(lexer*,fdm*,fieldint&);

    int column_pt1_count(lexer*,fdm*);
	int column_pt2_count(lexer*,fdm*);
	int column_pt3_count(lexer*,fdm*);
	int column_pt4_count(lexer*,fdm*);
    int column_pt4a_count(lexer*,fdm*);
    int column_pt6_count(lexer*,fdm*);


    void column_pt_resize(lexer*,fdm*);

	void sizeM_update(lexer*,fdm*);

    void fdm_update(fdm*);

// 2D CPT
    void cval2Dupdate1(lexer*,sliceint&);
    void cval2Dupdate2(lexer*,sliceint&);
    void cval2Dupdate4(lexer*,sliceint&);

    void column2D_pt1_update(lexer*,cpt2D&);
    void column2D_pt2_update(lexer*,cpt2D&);
    void column2D_pt4_update(lexer*,cpt2D&);

    void column2D_pt1(lexer*,sliceint&,cpt2D&);
    void column2D_pt2(lexer*,sliceint&,cpt2D&);
    void column2D_pt4(lexer*,sliceint&,cpt2D&);

    int column2D_pt1_count(lexer*);
    int column2D_pt2_count(lexer*);
    int column2D_pt4_count(lexer*);

    void cval_gcslpara1(lexer*,sliceint&);
    void cval_gcslpara2(lexer*,sliceint&);
    void cval_gcslpara4(lexer*,sliceint&);

    void cval2D_gcb1(lexer*, sliceint&);
    void cval2D_gcb2(lexer*, sliceint&);
    void cval2D_gcb4(lexer*, sliceint&);
    

    void sizeS_update(lexer*);

// solid update
	void solid_update(lexer*,fdm*);
	void gcsolid_gcb_remove(lexer*,fdm*);
	void gcsolid_gcb_seed(lexer*,fdm*);
	void gcsolid_gcb_dist(lexer*,fdm*);
	void gcsolid_buildflag(lexer*,fdm*, int&);
	void gcsolid_velflag1(lexer*,fdm*, int&);
	void gcsolid_velflag2(lexer*,fdm*, int&);
	void gcsolid_velflag3(lexer*,fdm*, int&);

// topo update
	void topo_update(lexer*,fdm*);
	void gcb_remove(lexer*,fdm*);
	void gcb_seed(lexer*,fdm*);
	void gcb_distbed(lexer*,fdm*);
	void gcb_buildflag(lexer*,fdm*, int**, int&);
	void gcb_velflag1(lexer*,fdm*, int **, int&);
	void gcb_velflag2(lexer*,fdm*, int **, int&);
	void gcb_velflag3(lexer*,fdm*, int **, int&);

	void velcell_update(lexer*, fdm*, int **, int ,double, double, double, int);
	void gctopo_pressureupdate(lexer*, fdm*, int **, int, field&);
	void gctopo_scalarupdate(lexer*, fdm*, int **, int, field&);



// 6DOF update
	void gcfb_update(lexer*,fdm*);
	void gcfb_buildflag(lexer*,fdm*, int**, int&);
	void gcfb_velflag1(lexer*,fdm*, int **, int&);
	void gcfb_velflag2(lexer*,fdm*, int **, int&);
	void gcfb_velflag3(lexer*,fdm*, int **, int&);
	void gcfb_seed(lexer*,fdm*);
	void gcfb_b_paraseed(lexer*,fdm*);
	void gcfb_x_paraseed(lexer*,fdm*);
	void gcfb_dist(lexer*,fdm*);
	void gcfb_velupdate(lexer*, fdm*, int **, int ,double, double, double, int);
	void gcfb_scalarupdate(lexer*, fdm*, int **, int, field&);
    void gcfb_update_extra_gcb(lexer*,fdm*,field&);
    void gcb_generic(lexer* p,field& f,int *gcb_count, int ***gcb);

// Forcing
    void forcing1(lexer*,fdm*,field&,field&,field&,field&,double);
	void forcing2(lexer*,fdm*,field&,field&,field&,field&,double);
	void forcing3(lexer*,fdm*,field&,field&,field&,field&,double);

// IBM
    void distguard(lexer*);
    void ccstateflag(lexer*,fdm*);
    void flagfield(lexer*);
    void flagfield_topo(lexer*);
    void tpflagfield(lexer*);
	void ndflag_update(lexer*);
    void flagbase(lexer*,fdm*);

// PARALLEL
    void gcparax(lexer*, field&, int);
    void gcparaxint(lexer*, fieldint&, int);
    void gcparax_test(lexer*, int);
    void gcparax_generic(lexer*, field&, int*, int***);
    void gcparacox_generic(lexer*, field&, int*, int***);
	void gcparaxvec(lexer*, vec&, int);
    void gcparax7(lexer*, double*, int);
    void gcparax7co(lexer*, double*, int);
	void gcparaxvec_sr(lexer*, vec&,cpt&,int);
	void gcparavoidx(lexer*, field&, int);
    void gcparax4a(lexer*, field&, int);
	void gcparacox(lexer*, field&, int);
	void verticalmax(lexer*,fdm*,double**);
    void verticalsum(lexer*,fdm*,double**);
    double timer();
    //Collective Communication
    void gather_int(int *, int, int *, int);
    void gather_double(double *, int, double *, int);
	void gatherv_int(int*, int, int*, int*, int*);
    void gatherv_double(double *, int, double *, int*, int*);
    void bcast_double(double *, int);
    //Utilities
    void walldistance(lexer*,fdm*,ghostcell*,convection*,reini*,ioflow*,field&);
	void walld_inflow(lexer*,fdm*,ghostcell*,double*);
	void walld_outflow(lexer*,fdm*,ghostcell*,double*);
    void gcwait(lexer*);
    void gcwait7(lexer*);

    MPI_Comm mpi_comm;

// Slice
    // epol
    void gcsl_start1(lexer*,slice&, int);
	void gcsl_start2(lexer*,slice&, int);
	void gcsl_start3(lexer*,slice&, int);
	void gcsl_start4(lexer*,slice&, int);
	void gcsl_start4a(lexer*,slice&, int);

    void gcsl_start4int(lexer*,sliceint&, int);

    void gcsl_start4V(lexer*,vec2D&, int, cpt2D&);
    
    void gcsl_solidupdate(lexer*);

    void gcsldistro1(lexer*, slice&,int, int, int, double, int, int, int);
	void gcsldistro2(lexer*, slice&,int, int, int, double, int, int, int);
	void gcsldistro3(lexer*, slice&,int, int, int, double, int, int, int);
	void gcsldistro4(lexer*, slice&,int, int, int, double, int, int, int);
	void gcsldistro4a(lexer*, slice&,int, int, int, double, int, int, int);
    
    void gcsldistro4V(lexer*, vec2D&, cpt2D &, int, int, int, double, int, int, int,int);

    void gcsldistro4int(lexer*, sliceint&,int, int, int, double, int, int, int);

    int gcsleval1(lexer*,int,int,int);
	int gcsleval2(lexer*,int,int,int);
	int gcsleval3(lexer*,int,int,int);
	int gcsleval4(lexer*,int,int,int);
	int gcsleval4a(lexer*,int,int,int);

	void gcsl_tpflag(lexer*);

    double gcsl_ipol1(lexer*,slice&);
    double gcsl_ipol2(lexer*,slice&);
	double gcsl_ipol1a(lexer*,slice&);
    double gcsl_ipol2a(lexer*,slice&);
    double gcsl_ipol4(lexer*,slice&);
    double gcsl_ipol4eta(lexer*,slice&,slice&);
    double gcsl_ipolint(lexer*,sliceint&);

	double gcsl_ccipol4(lexer*,slice&,double,double);
	double gcsl_lint4(lexer*,slice&,int&,int&,double,double);

    void gcsl_setbc1(lexer*);
    void gcsl_setbc2(lexer*);
    void gcsl_setbc4(lexer*);
    void gcsl_setbcio(lexer*);

    // Slice BCs
    void gcsl_neumann(slice&,int,int,int);
    void gcsl_neumann_x(slice&,int,int,int);
    void gcsl_neumann_int(sliceint&,int,int,int);
	void gcsl_noslip(slice&,int,int,int);
    void gcsl_sommerfeld(lexer*,slice&,int,int,int);
    void gcsl_outflow(lexer*,slice&,int,int,int);
    void gcsl_outflow_fsf(lexer*,slice&,int,int,int);
    
    void gcsl_neumannV(vec2D&,cpt2D&,int,int,int,int);
    void gcsl_neumannV_x(vec2D&,cpt2D&,int,int,int,int);
	void gcsl_noslipV(vec2D&,cpt2D&,int,int,int,int);


    // parallel
    void gcslparax(lexer*, slice&, int);
    void gcslparax_fh(lexer*, slice&, int);
    void gcslparax_int(lexer*, sliceint&, int);
    void gcslparacox(lexer*, slice&, int);
    void gcslparacox_int(lexer*, sliceint&, int);
    void gcslwait(lexer*);
    void gcslflagx(lexer*, int*);
    void gcparaxvec2D(lexer*, vec2D&, int,cpt2D&);
    void gcslparaxvec_slr(lexer*, vec2D&,cpt2D&,int);
    void gcxslupdate(lexer*);


    double mini1(fdm*,lexer*, field&);
    double maxi1(fdm*,lexer*, field&);

    double mini2(fdm*,lexer*, field&);
    double maxi2(fdm*,lexer*, field&);

    double mini3(fdm*,lexer*, field&);
    double maxi3(fdm*,lexer*, field&);

    double mini4(fdm*,lexer*, field&);
    double maxi4(fdm*,lexer*, field&);

    int ii,jj,kk;
    int ic,jc,kc;

    double dx;

    int imin,imax,jmax,jmin,kmin,kmax;



	void gcdistro1(lexer *p,field&,int, int, int, int, double, int, int, int);
	void gcdistro2(lexer *p,field&,int, int, int, int, double, int, int, int);
	void gcdistro3(lexer *p,field&,int, int, int, int, double, int, int, int);
	void gcdistro4(lexer *p,field&,int, int, int, int, double, int, int, int);
	void gcdistro4a(lexer *p,field&,int, int, int, int, double, int, int, int);

	void gcdistro4V(lexer *p,fdm*, vec&, int, int, int, double, int, int, int, int);
    void gcdistro4aV(lexer *p,fdm*, vec&, int, int, int, double, int, int, int, int);

	int gceval1(lexer*,int,int,int);
	int gceval2(lexer*,int,int,int);
	int gceval3(lexer*,int,int,int);
	int gceval4(lexer*,int,int,int);
	int gceval4a(lexer*,int,int,int);

    void nse1(lexer*, fdm*, field&, int);
    void nse2(lexer*, fdm*, field&, int);
    void nse3(lexer*, fdm*, field&, int);
    void nse4(lexer*, fdm*, field&, int);

    void nse1_conv(lexer*, fdm*, field&, int, double);
    void nse2_conv(lexer*, fdm*, field&, int, double);
    void nse3_conv(lexer*, fdm*, field&, int, double);

	void dirichlet_para(field&,double,int,int,int);
	void dirichlet_ortho(lexer*,field&,double,int,int,int);
	void neumann(field&,int,int,int);
    void gcb_debug(field&,int,int,int);
	void neumann_press(lexer*,field&,double,int,int,int);
	void extend(field&,double,int,int,int);
	void largeval(field&,double,int,int,int);
	void largevaladd(field&,double,int,int,int);
	void outflow(lexer*,field&,int,int,int);
    void sommerfeld(lexer*,field&,int,int,int);
	void inflowbc(field&,double,int,int,int);
    void potentialbc(lexer*,field&,int,int);
    void neumann_all(field&,int,int,int);
    void extend_all(field&,double,int,int,int);
    void lsm(lexer*,field&,double,int,int,int);
    void noslip(field&,double,int,int,int);
    void imagepoint(field&, double&, double&,double,int);
	void atmosphere(lexer*,field&,int,int,int);
    void heatbc(lexer*,field&,int,int,int);
	void fbvel1(lexer*,field&,double,int,int,int);
	void fbvel2(lexer*,field&,double,int,int,int);
	void fbvel3(lexer*,field&,double,int,int,int);
	void fbpress(lexer*,field&,double,int,int,int);
	void gravity_press(lexer*,field&,double,int,int,int);
    void kinematic_bed(lexer*,field&,double,int,int,int);
    void fivec(lexer*,double*,sliceint&);
    void fivec2D(lexer*,double*,sliceint&);
    void fivec_vel(lexer*,double*,sliceint&);
    void fivec2D_vel(lexer*,double*,sliceint&);

	void gcV_neumann(vec&,int,int,int,int);
	void gcV_lsm(lexer*,vec&, double,int,int,int,int);
    void gcV_neumann_all(vec&, int,int,int,int);


private:
    const int size;
    const int tag1,tag2,tag3,tag4,tag5,tag6;
	double deltax;
	int margin, paramargin;
	double  y[15],dP[15], x[15],pos[15];
	double val[10];
	int m,q,qq,qn,g;
	int bc_label;
	double wallvalue,x_ip,val_ip,gamma;
	int orderdir,orderdir2,orderext,orderext2,orderpress;
	double Qi,weight;
	int close;
	double dist;
    int count,check;
    double starttime,endtime;
    const double eps;
    int offset,ys;
    int gcval_topodist;
	int gclabel_outflow;
    int gclabel_u, gclabel_v, gclabel_w, gclabel_k, gclabel_e;
    int gclabel_utopo, gclabel_vtopo, gclabel_wtopo;
    int gclabel_u_orth,gclabel_v_orth,gclabel_w_orth,gclabel_press,gclabel_lsm;
    int gclabel_u_in,gclabel_v_in,gclabel_w_in,gclabel_press_in,gclabel_lsm_in;
	int gc_pressure_extend;
	int gclabel_u_out, gclabel_v_out, gclabel_w_out;
	int gclabel_vel;
	int rank;
	int nb[6],stag[6],rtag[6];
	int **isend,**irecv;
	double **dsend,**drecv;
	double *trecv;

	density_f *pdens;

    double originx,originy,originz;

// PARALLEL


	MPI_Request sreq1,sreq2,sreq3,sreq4,sreq5,sreq6;
	MPI_Request rreq1,rreq2,rreq3,rreq4,rreq5,rreq6;

	MPI_Request sreq[6],rreq[6];


	MPI_Status status;


	int tag;
	double **send,**recv;
	double *send1,*send2,*send3,*send4,*send5,*send6;
	double *recv1,*recv2,*recv3,*recv4,*recv5,*recv6;
	int *isend1,*isend2,*isend3,*isend4,*isend5,*isend6;
	int *irecv1,*irecv2,*irecv3,*irecv4,*irecv5,*irecv6;
	int precv[6];
	double recvsum,recvmin,recvmax;
	int recvisum,recvimin,recvimax;
	int awa_label;
	const int gcx;
	int gcx_count[6];

	double epphi;

    double v1,v2,v3,v4;
    double wa,wb;
    double x1,x2;
    double value;

// 6DOF
	int ***gcbfb,*gcbfb_count;
	int ***gcxfb,*gcxfb_count;

    fdm *a;

};
#endif
