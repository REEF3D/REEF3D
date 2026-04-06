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

#ifndef LEXER_H_
#define LEXER_H_

#include<iostream>
#include<cstdlib>
#include<iomanip>
#include<math.h>
#include"resize.h"
#include"increment.h"
#include"position.h"
#include"interpolation.h"
#include<fstream>
#include"looping.h"
#include<vector>
#include"grid.h"
#include"control.h"
#include"coordinates.h"

class weno_nug_func;
class ghostcell;

using namespace std;

class lexer : virtual public resize_class, public grid, public control,
              public position, public interpolation, public coordinates
{
public:

	lexer();
	virtual ~lexer();

//-----functions------------------
//---- setup
    void lexer_read(ghostcell*);
    void flagini();
	void gridini(ghostcell*);
    void gcd_ini(ghostcell*);
    void makeflag(int*);

	void read_grid();
	void control_calc();
	void gridsize();
	void vecsize(ghostcell*);
    void gcbextra_est(ghostcell*);
	void vellast();
	void lexer_ini();
	int conv(double);

    // 2D
    void grid2Dsize();
    void flagini2D();
	void gridini2D();


//-----data-----------------------

	//REEF3D

	int pointnum,cellnum,tpcellnum;
	int cellnum1,cellnum2,cellnum3;
    int pointnumtot,cellnumtot;
    int N4,N4_row,N4_col;
    int N7,N7_row,N7_col;
	int surf_tot;
	int *flag1,*flag2,*flag3,*flag4,*flag5,*flag7,*flag;
    int *flagsf1,*flagsf2,*flagsf3,*flagsf4;

    // flag
    double *flag_solid,*flag_topo;
    double *data;
	double *topobed,*solidbed,*bed,*depth;
    int *wet,*wet_n;
    int *deep;
    int gcbextra;
    int solidread,toporead,porousread,topoforcing;


    //GHOSTCELL
	int **gcb1,**gcb2,**gcb3,**gcb4,**gcb4a;
	int **gcin, **gcout, **gcpress;
	int **gcin4a, **gcout4a;
	double *gcd1,*gcd2,*gcd3,*gcd4,*gcd4a;
	double **gcn;
	int gcextra1,gcextra2,gcextra3,gcextra4,gcextra4a,gcextra6;

    int gcdf1_count,gcdf2_count,gcdf3_count,gcdf4_count;
    int **gcdf1,**gcdf2,**gcdf3,**gcdf4;
    int gcsldfeta4_count,gcsldfbed4_count;
    int **gcsldfeta4,**gcsldfbed4;

	int gcwall_count, gcin_count, gcout_count, gcpress_count, gcfsf_count, gcbed_count;
	int gcin4a_count, gcout4a_count;
	int gcb1_count,gcb2_count,gcb3_count,gcb4_count,gcb4a_count;
	int gcpara_sum, gcparaco_sum;
	int gcb_fix,gcb_solid,gcb_topo,gcb_fb, solid_gcb_est, topo_gcb_est, solid_gcbextra_est, topo_gcbextra_est, tot_gcbextra_est;
	int gcb_sediment_est, gcb_floating_est;
    int bcside1,bcside2,bcside3,bcside4,bcside5,bcside6;

    // serial periodic BC
    int periodic1,periodic2,periodic3;
    int periodicX1,periodicX2,periodicX3,periodicX4,periodicX5,periodicX6;

    int **dgc1,**dgc2,**dgc3,**dgc4;
    int dgc1_count,dgc2_count,dgc3_count,dgc4_count;

	// PARALLEL
	int** gcpara1;
	int** gcpara2;
	int** gcpara3;
	int** gcpara4;
	int** gcpara5;
	int** gcpara6;

	int** gcparaco1;
	int** gcparaco2;
	int** gcparaco3;
	int** gcparaco4;
	int** gcparaco5;
	int** gcparaco6;

    int*** gcx7;
    int* gcx7_count;
    int*** gcxco7;
    int* gcxco7_count;

	int gcpara1_count, gcpara2_count, gcpara3_count, gcpara4_count, gcpara5_count, gcpara6_count;
	int gcparaco1_count, gcparaco2_count, gcparaco3_count, gcparaco4_count, gcparaco5_count, gcparaco6_count;
    int gcslpara1_count, gcslpara2_count, gcslpara3_count, gcslpara4_count;
    int gcslparaco1_count, gcslparaco2_count, gcslparaco3_count, gcslparaco4_count;
	int nb1,nb2,nb3,nb4,nb5,nb6;
	int mx,my,mz;
	int mpi_edgenum,mpi_nodes,mpi_size;
	int *mpi_index, *mpi_edges;

	int ulast,vlast,wlast,flast,ulastsflow;
	int velcorr;
	int stencil;

	// Solver
    int *range_col4,*range_row4,*range_col7,*range_row7;
	int *sizeM1,*sizeM2,*sizeM3,*sizeM4,*sizeM4a,*sizeM6,*sizeM9;
    int *sizeS1,*sizeS2,*sizeS4;
	int mglevel_max,*MGL;

	// SMO
	int veclength;
    int C4_size,C4a_size,C6_size;
    int C1_2D_size,C2_2D_size,C4_2D_size;
    int M_size,M_2D_size;

    //SLICE
    int *flagslice1,*flagslice2,*flagslice4;
    int *mgcsl1,*mgcsl2,*mgcsl3,*mgcsl4,*mgcsl4a;
    int ***gcslorig1,***gcslorig2,***gcslorig3,***gcslorig4,***gcslorig4a;
	int gcsldirsize1,gcsldirsize2,gcsldirsize3,gcsldirsize4,gcsldirsize4a;

    int slicenum,vec2Dlength;

    int pointnum2D,cellnum2D,cellnumtot2D,polygon_sum;

    // SLICE ghostcell
    int gcbsl1_count,gcbsl2_count,gcbsl3_count,gcbsl4_count,gcbsl4a_count;
    int gcslin_count,gcslout_count;
    int gcslawa1_count,gcslawa2_count;
    int **gcbsl1,**gcbsl2,**gcbsl3,**gcbsl4,**gcbsl4a;
	int **gcslin, **gcslout;
    int **gcslawa1, **gcslawa2;

    int gcsl_extra1,gcsl_extra2,gcsl_extra3,gcsl_extra4,gcsl_extra4a;

	int **dgcsl1,**dgcsl2,**dgcsl3,**dgcsl4;
	int dgcsl1_count,dgcsl2_count,dgcsl3_count,dgcsl4_count;

    int **ggcsl1,**ggcsl2,**ggcsl3,**ggcsl4,**ggcsl4a;
    int *ggcslmem1,*ggcslmem2,*ggcslmem3,*ggcslmem4,*ggcslmem4a;
    int ggcslcount1,ggcslcount2,ggcslcount3,ggcslcount4,ggcslcount4a;
    int ggcslsize1,ggcslsize2,ggcslsize3,ggcslsize4,ggcslsize4a;

    // SLICE parallel
	int** gcslpara1;
	int** gcslpara2;
	int** gcslpara3;
	int** gcslpara4;

	int** gcslparaco1;
	int** gcslparaco2;
	int** gcslparaco3;
	int** gcslparaco4;


    // flow parameters
    const double cmu;
    double deltax,sigT,Ui,Ua,Uo;
    double Ho,Hi;

    // 6DOF
	double ufb,vfb,wfb;
	double pfb,qfb,rfb;
	double ufbi,vfbi,wfbi;
	double pfbi,qfbi,rfbi;
	double xg,yg,zg;
	double xgn,ygn,zgn;
	double phi_fb,theta_fb,psi_fb;
	double ufbmax, vfbmax, wfbmax;
    int mooring_count, net_count;

    // FSI
    int FSI_count;

	// time + iterations
	int inneriter,count,solveriter,preconiter,count_statestart;
    int solver_status,solver_error;
    int sediter;
    double final_res;
	double dt,dt_old,simtime,viscmax;
	double mindt,maxdt;
	double umax,vmax,wmax,epsmax,kinmax,pressmin,pressmax,omegamax;
	double presstime,veltime,reinitime,turbtime,plstime,itertime;
	double sedsimtime,sedwavetime;
	double wavecalctime;
	double meantime,totaltime;
	double gcmeantime,gctotaltime;
	double Xmeantime,Xtotaltime;
	double maxbed, minbed;
	double susptime,maxtopovel;
	double gctime, xtime;
	double volume1,volume2,volume3;
	double Qi,Qo;
	double dtsed,sedtime,slidecells;
	double bedmax,bedmin;
	double field4time;
    double printtime, sedprinttime,fsfprinttime,fsfsedprinttime,probeprinttime,stateprinttime,exportprinttime;
    double wavetime;
    double RK_alpha;

	// solver watch
	int uiter,viter,witer;
	int kiniter,epsiter;
	int poissoniter, laplaceiter;
	int lsmiter;
	int suspiter,topoiter;
	int heatiter,concentrationiter;
	int printcount, printcount_sixdof;
	double utime,vtime,wtime;
    double recontime,fsftime;
    double dftime;
	double kintime,epstime;
	double poissontime, laplacetime, matrixtime, ptime;
    double sftime,fbtime,fsitime;
    double fbdt,fbmax;
    double sfdt,sfmax;
	double lsmtime,heattime,concentrationtime;
	double printouttime;
	double phimean,phiout,phiin;
    double fsfin,fsfout;
    double fsfinval,fsfoutval;
	double pcnorm,ucnorm,vcnorm,wcnorm;
    double alpha;
    double pressgage;

	// wave coefficients
	double wT,wV,wH,wA,wL,wd,ww,wk,wC;
	double wHs,wAs,wwp,ww_s,ww_e,wTp,wLp;
	int wN;
    double wts,wte;

    // free surface
    double psi,psi0;
	int pressval;

// PARALELL
    int mpirank;
	int gcx_1range1[7],gcx_3range1[7];
	int gcx_1range2[7],gcx_3range2[7];
	int gcx_1range3[7],gcx_3range3[7];
	int gcx_1range4[7],gcx_3range4[7];

    weno_nug_func *wenofunc;

// sigma coordinate
    double *sig;
    double *sigx,*sigy,*sigz,*sigt;
    double *sigxx;
};

#endif
