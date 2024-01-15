/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

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

#ifndef LEXER_H_
#define LEXER_H_

class weno_nug_func;
class ghostcell;

using namespace std;

class lexer : public increment, public resize_class, public position, public interpolation
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
    void gridini_patchBC();
    void makeflag(int*);
    void sigma_coord_ini();
	
	void read_grid();
	void read_control();
	void control_calc();
	void ini_default();
	void assign_margin();
	void ctrlsend();
	void ctrlrecv();
	int maxparacount();
	void gridsize();
	void vecsize(ghostcell*);
    void gcbextra_est(ghostcell*);
	void vellast();
	void indices_minmax();
	void lexer_ini();
    void lexer_gridspacing(ghostcell*);
	void parse();
	void fieldlogic();
	int conv(double);
    
    // 2D
    void grid2Dsize();
    void flagini2D();
	void gridini2D();


//-----data-----------------------

	//REEF3D
	double dx,dy,dz;
    double *xpoint,*ypoint,*zpoint;
    double *xnode,*ynode,*znode;
    
    
	int imin,imax,jmin,jmax,kmin,kmax;
    int kmaxF;
	int pointnum,cellnum,tpcellnum;
	int cellnum1,cellnum2,cellnum3;
    int pointnumtot,cellnumtot;
    int N4,N4_row,N4_col;
    int N7,N7_row,N7_col;
	double originx,originy,originz;
    double endx,endy,endz;
	double global_xmin,global_ymin,global_zmin;
	double global_xmax,global_ymax,global_zmax;
	int origin_i, origin_j, origin_k;
	int gknox,gknoy,gknoz;
	int surf_tot;
	int *flag1,*flag2,*flag3,*flag4,*flag5,*flag7,*flag;
    int *flagsf1,*flagsf2,*flagsf3,*flagsf4;
    int *BC;
	int*mgflag;
    double *flag_solid,*flag_topo;
    double *data;
	double *topobed,*solidbed,*bed,*depth;
    int *wet,*wet_n;
    int *deep;
	int *tpflag,*ndbaseflag;
	int *mgc1,*mgc2,*mgc3,*mgc4,*mgc4a,*mgc6;
	int ***gcorig1,***gcorig2,***gcorig3,***gcorig4,***gcorig4a,***gcorig6;
	int gcdirsize1,gcdirsize2,gcdirsize3,gcdirsize4,gcdirsize4a,gcdirsize6;
	int i_dir,j_dir,k_dir;
	double x_dir,y_dir,z_dir;
    int gcbextra;
    int solidread,toporead,porousread;


    //GHOSTCELL
	int **gcb1,**gcb2,**gcb3,**gcb4,**gcb4a,*gcb6;
	int **gcin, **gcout, **gcpress,**gcin6, **gcout6;
	int **gcin4a, **gcout4a;
	double *gcd1,*gcd2,*gcd3,*gcd4,*gcd4a;
	double **gcn;
	int *gcside4;
	int gcside4_size;
	int gcextra1,gcextra2,gcextra3,gcextra4,gcextra4a,gcextra6;
    
    int gcdf1_count,gcdf2_count,gcdf3_count,gcdf4_count;
    int **gcdf1,**gcdf2,**gcdf3,**gcdf4;

	int **dgc1,**dgc2,**dgc3,**dgc4;
	int dgc1_count,dgc2_count,dgc3_count,dgc4_count;

	int gcwall_count, gcin_count, gcout_count, gcpress_count, gcfsf_count, gcbed_count;
    int gcin6_count, gcout6_count;
	int gcin4a_count, gcout4a_count;
	int gcb1_count,gcb2_count,gcb3_count,gcb4_count,gcb4a_count;
	int gcpara_sum, gcparaco_sum;
	int gcb_fix,gcb_solid,gcb_topo,gcb_fb, solid_gcb_est, topo_gcb_est, solid_gcbextra_est, topo_gcbextra_est, tot_gcbextra_est;
	int gcb_sediment_est, gcb_floating_est;
    int bcside1,bcside2,bcside3,bcside4,bcside5,bcside6;
    
    // serial periodic BC
    int periodic1,periodic2,periodic3;
    int periodicX1,periodicX2,periodicX3,periodicX4,periodicX5,periodicX6;
    
    int **gc4periodic;
    int **gc4aperiodic;
    int *gc4periodic_count;
    int *gc4aperiodic_count;
    int gc4periodic_maxcount;
    
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
	int maxpara;
	int nb1,nb2,nb3,nb4,nb5,nb6;
	int mx,my,mz;
    int mi,mj,mk;
	int mpi_edgenum,mpi_nodes,mpi_size;
	int *mpi_index, *mpi_edges;
	
	int ulast,vlast,wlast,flast,ulastsflow;
	int velcorr;
	int* ictrl;
	double* dctrl;
	int ctrlsize;
	int stencil;	

	// Solver
	int *colnum;
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
    int *flagslice1,*flagslice2,*flagslice4,*tpflagslice;
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
	double *gcdsl1,*gcdsl2,*gcdsl3,*gcdsl4,*gcdsl4a;


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

    // Hydrodynamics Models
    int A10;
    
    // SFLOW
	int A209,A210,A211,A212,A214,A215,A216,A217,A218,A219,A220,A221,A230,A240,A241,A242,A243,A244,A245,A246,A248;
    int A251,A260;
    double A261,A262;
    double A223,A244_val,A245_val,A247,A249,A251_val;
    double A250;
    
    // FNPF
    int A310,A311,A312,A313,A320,A321,A322,A323,A329,A343,A344,A345,A347,A348;
    double A340,A341,A342,A344_val,A345_val,A346;
    int A350,A351,A352,A353,A357,A361,A362,A363,A368;
    double A354,A355,A356,A365; 
    
    // NSEWAVE
    int A410;
    double A440;
    
    // NHFLOW
    int A501,A510,A511,A512,A514,A515,A516,A517,A518;
    int A520,A521;
    double A522,A523;
    double A531;
    int A540,A543;
    double A541,A542,A544;
    int A550,A551,A552,A553;
    
	// boundary conditions
	int B10,B20,B23;
    int B30,B32,B33;
    double B31,B32_x,B32_y,B32_z;    
    int B60,B61,B71,B75,B76,B77,B84,B85,B81,B82,B86,B87,B89,B90,B91,B92,B93,B94,B98,B99,B101,B105,B106,B107;
	int B136,B138,B138_1,B138_2,B139;
    int B180,B191,B192,B240,B241,B242,B243;
	double B29,B50,B51,B52,B53,B54,B55,B56,B81_1,B81_2,B81_3,B83,B117,B87_1,B87_2,B88;
	double B91_1,B91_2,B93_1,B93_2,B94_wdt,B96_1,B96_2,B102,B105_1,B105_2,B105_3;
	double *B71_val,*B71_dist,*B71_b,*B71_x,*B71_y;
	double *B106_b,*B106_x,*B106_y;
    double *B107_xs,*B107_xe,*B107_ys, *B107_ye, *B107_d;
    int B108;
    double *B108_xs,*B108_xe,*B108_ys, *B108_ye, *B108_d;
    int B110;
    double B110_zs,B110_ze;
	double B111_zs,B111_ze;
    double B112_zs,B112_z2,B112_ze;
    int B115,B116,B125,B127;
    double B120,B122,B123,B125_y;
    int B130,B133;
    double B131,B132_s,B132_e;
    double B134,B135;
    int B160, B170;
    int B181,B182,B183;
	double B181_1,B181_2,B181_3,B182_1,B182_2,B182_3,B183_1,B183_2,B183_3;
	double B191_1,B191_2,B191_3,B191_4,B192_1,B192_2,B192_3,B192_4;
	double B194_s,B194_e;
    
    int B411,B412,B413,B414,B415,B416,B417,B418,B421,B422;
    int *B411_ID;
    double *B411_Q;
    int *B412_ID;
    double *B412_pressBC;
    int *B413_ID;
    double *B413_h;
    int *B414_ID;
    double *B414_Uio;
    int *B415_ID;
    double *B415_U,*B415_V,*B415_W;
    int *B416_ID;
    double *B416_alpha;
    int *B417_ID;
    double *B417_Nx,*B417_Ny,*B417_Nz;
    int *B418_ID;
    int *B418_pio;
    int *B421_ID,*B421_Q;
    int *B422_ID,*B422_FSF;
    int B440;
    int *B440_ID,*B440_face;
    double *B440_xs,*B440_xe,*B440_ys,*B440_ye;
    int B441;
    int *B441_ID,*B441_face;
    double *B441_xs,*B441_xe,*B441_ys,*B441_ye,*B441_zs,*B441_ze;
    int B442;
    int *B442_ID,*B442_face;
    double *B442_xm,*B442_ym,*B442_zm,*B442_r;
    
	double *B240_D, *B240_C, *B240_xs, *B240_xe, *B240_ys, *B240_ye, *B240_zs, *B240_ze;
    double B260,B264,B267;
    int B269,B270;
    double *B270_xs, *B270_xe, *B270_ys, *B270_ye, *B270_zs, *B270_ze, *B270_n, *B270_d50, *B270_alpha, *B270_beta;
    int B274;
    double *B274_xc,*B274_yc,*B274_zs,*B274_ze,*B274_r, *B274_n, *B274_d50, *B274_alpha, *B274_beta;
    int B281;
    double *B281_xs, *B281_xe, *B281_ys, *B281_ye, *B281_zs, *B281_ze, *B281_n, *B281_d50, *B281_alpha, *B281_beta;
    int B282;
    double *B282_xs, *B282_xe, *B282_ys, *B282_ye, *B282_zs, *B282_ze, *B282_n, *B282_d50, *B282_alpha, *B282_beta;
	int B291;
    double *B291_xs, *B291_xe, *B291_ys, *B291_ye, *B291_zs, *B291_ze, *B291_d, *B291_n, *B291_d50, *B291_alpha, *B291_beta;
    int B295;
    int B308,B310,B311;
    double B309;
    double *B310_xs, *B310_xe, *B310_ys, *B310_ye, *B310_zs, *B310_ze, *B310_N, *B310_D, *B310_Cd;
    double *B311_xm, *B311_ym, *B311_r, *B311_zs, *B311ze, *B311_N, *B311_D, *B311_Cd;
    int B321;
    double *B321_xs, *B321_xe, *B321_ys, *B321_ye, *B321_zs, *B321_ze, *B321_N, *B321_D, *B321_Cd;
    int B322;
    double *B322_xs, *B322_xe, *B322_ys, *B322_ye, *B322_zs, *B322_ze, *B322_N, *B322_D, *B322_Cd;
	
    // Concentration Options
	double C1,C2,C3,C4,C5;
	int C9,C10,C15,C20;
	double C50_1,C50_2;
	double C51,C52,C53,C54,C55,C56;
	double C57_1,C57_2,C57_3,C57_4;
	double C58_1,C58_2,C58_3,C58_4;
	int C75;
	double *C75_x,*C75_z,*C75_a,*C75_s,*C75_l,*C75_v;

	// discretization
	int D10,D11,D20,D21,D30,D31,D33,D37;

	// Free Surface
	int F10,F30,F31,F32,F34,F35,F36,F40,F44,F46,F47,F49,F50,F150,F151;
	double F33,F39,F42,F43,F45;
	double F51,F52,F53,F54,F55,F56;
    int F50_flag;
	double F57_1,F57_2,F57_3,F57_4;
	double F58_1,F58_2,F58_3,F58_4;
    double F59_xm, F59_ym, F59_zs, F59_ze, F59_r;
	double F60,F61,F62,F63;
	int F64;
	double F64_xs,F64_ys,F64_zs,F64_alpha;
	int F70;
	double *F70_xs, *F70_xe, *F70_ys, *F70_ye, *F70_zs, *F70_ze;
	int F71;
	double *F71_xs, *F71_xe, *F71_ys, *F71_ye, *F71_zs, *F71_ze;
	int F72;
	double *F72_xs, *F72_xe, *F72_ys, *F72_ye, *F72_h;
	int F80,F85;
	double F84;
    
    int F300,F305,F310,F350;
	double F321,F322,F323,F360,F361,F362;
	int F369,F370,F371,F374,F375,F378,F379;
    double *F369_x,*F369_z,*F369_a,*F369_s,*F369_l,*F369_v;
	double *F370_xs, *F370_xe, *F370_ys, *F370_ye, *F370_zs, *F370_ze;
	double *F371_xs, *F371_xe, *F371_ys, *F371_ye, *F371_zs, *F371_ze;
	double *F374_xc, *F374_zc, *F374_r;
    double *F375_xc, *F375_zc, *F375_r;
    double *F378_xc, *F378_yc,*F378_zc, *F378_r;
    double *F379_xc, *F379_yc,*F379_zc, *F379_r;
	double F380,F381,F382;
	int F390,F391,F394,F395,F398,F399;
	double *F390_xs, *F390_xe, *F390_ys, *F390_ye, *F390_zs, *F390_ze;
	double *F391_xs, *F391_xe, *F391_ys, *F391_ye, *F391_zs, *F391_ze;
    double *F394_xc, *F394_zc, *F394_r;
    double *F395_xc, *F395_zc, *F395_r;
    double *F398_xc, *F398_yc,*F398_zc, *F398_r;
    double *F399_xc, *F399_yc,*F399_zc, *F399_r;
    
	// Grid Options
    int G1,G2,G3;
	int G10,G11,G12,G20,G21,G22,G30;
	int G40;

	// Heat Options
	double H1,H2;
	int H3,H4,H9,H10,H15;
	double H4_beta1,H4_beta2,H50_1,H50_2;
	double H51,H52,H53,H54,H55,H56;
	double H57_1,H57_2,H57_3,H57_4;
	double H58_1,H58_2,H58_3,H58_4;
    int H61,H62,H63,H64,H65,H66;
    double H61_T,H62_T,H63_T,H64_T,H65_T,H66_T;
	
	// Initialize Options
	int I10,I11,I12,I13,I30,I40,I41,I44,I56;
	double I21,I50,I55,I58_1,I58_2;
    int I230;
    double I231,I232,I233;
    int I240;
    double I241;

	// Numerical Options
	int N10,N11,N40,N45,N46,N48,N60;
	double N41,N43,N44,N47,N49,N50,N61;

	// MPI Options
	int M10;

	// Print options
	int P10,P11,P12,P14,P15,P18,P20,P21,P23,P24,P25,P26,P27,P28,P29,P35,P40,P41,P43,P44,P45,P50,P51,P52,P53,P54,P56,P57,P58,P59;
	int P61,P62,P63,P64,P66,P67,P68,P71,P72,P73,P74,P75,P76,P77,P78,P79,P81,P82,P85,P92,P101,P120,P121,P122,P123,P124,P125,P126;
	int P150,P151,P152,P180,P181,P184,P185,P190,P191,P194,P195,P351,P352;
	double P22,P30,P34,P42;
	double *P35_ts,*P35_te,*P35_dt;
    double P43_xs,P43_xe,P43_ys,P43_ye;
	double *P50_x,*P50_y;
	double *P51_x,*P51_y;
	double *P52_y,*P56_x;
	double P55;
    double *P58_x,*P58_y,*P58_T;   
	double *P61_x,*P61_y,*P61_z;
	double *P62_xs,*P62_ys,*P62_zs,*P62_xe,*P62_ye,*P62_ze;
    double *P63_x,*P63_y;
    double *P64_x,*P64_y,*P64_z;
	double *P67_x;
    double *P68_x,*P68_zs,*P68_ze;
	double *P81_xs,*P81_xe,*P81_ys,*P81_ye,*P81_zs,*P81_ze;
	double *P85_x,*P85_y,*P85_r,*P85_cd,*P85_cm;
	double P91;
	double P101_xm,P101_ym,P101_zs,P101_ze,P101_r1,P101_r2;
    int P110;
    double P111;
	double *P121_x,*P121_y;
	double *P123_y,*P124_x;
	double *P125_x,*P125_y;
	double P182;
    int *P184_its,*P184_ite,*P184_dit;
    double *P185_ts,*P185_te,*P185_dt;
    double P192;
    int *P194_its,*P194_ite,*P194_dit;
    double *P195_ts,*P195_te,*P195_dt;
    int P230,P240;
    double *P230_x,*P240_x;
	double *P351_x,*P351_y;
	double *P352_x,*P352_y;
    
    // Particles
    int Q10,Q24,Q29,Q43;
    double Q21,Q22,Q23,Q25;
    double Q31;
    double Q41;
    
    int Q101,Q110;
    double Q102;
    double *Q110_xs,*Q110_xe,*Q110_ys,*Q110_ye,*Q110_zs,*Q110_ze;
    int Q111,Q112,Q113;
    double Q111_x,Q112_y,Q113_z;
    
    int Q180,Q181;
    double Q182;
    

	// Sediment Transport
	int S10,S11,S12,S15,S16,S17,S23,S27,S32,S33,S34,S37,S41,S42,S43,S44,S50,S60,S73,S77,S78,S79,S80,S83,S84,S90,S91,S100,S101;
	double S13,S14,S19,S20,S21,S22,S23_val,S24,S26_a,S26_b,S30,S45,S46,S47,S48,S57,S71,S72,S81,S82,S92,S93,S116;
	double *S73_val,*S73_dist,*S73_b,*S73_x,*S73_y;
    double S77_xs,S77_xe;

	// Turbulence
	int T10,T11,T12,T21,T33,T36,T39,T41;
	double T31,T32,T35,T37,T38,T42,T43;

	// Waterflow
	double W1,W2,W3,W4,W5,W6,W7,W10,W_fb;
    int W11,W12,W13,W14,W15,W16;
    double W11_u,W11_v,W11_w,W12_u,W12_v,W12_w,W13_u,W13_v,W13_w,W14_u,W14_v,W14_w,W15_u,W15_v,W15_w,W16_u,W16_v,W16_w;
    double W20,W21,W22,W31,W29_x,W29_y,W29_z;
	int W30;
    int W41;
    double *W41_xc,*W41_yc,*W41_zs,*W41_ze,*W41_vel,*W41_beta;
    double W50;
    int W50_air;
    int W90;
    double W95,W96,W97,W98;
    int W101;
    double W102_c,W102_phi;
    double W103,W104;
    int W110,W111;
    double W112;
    
    // 6DOF
	double ufb,vfb,wfb;
	double pfb,qfb,rfb;
	double ufbi,vfbi,wfbi;
	double pfbi,qfbi,rfbi;
	double ufbn,vfbn,wfbn;
	double pfbn,qfbn,rfbn;
	double xg,yg,zg;
	double xgn,ygn,zgn;
	double phi_fb,theta_fb,psi_fb;
	double ufbmax, vfbmax, wfbmax;
	//Eigen::Matrix3d quatRotMat;	
    int X10,X12,X14,X15,X19,X11_u,X11_v,X11_w,X11_p,X11_q,X11_r,X21,X22,X23,X24,X31,X32,X33,X34,X38;
    int X39,X40,X45,X46,X47,X48,X49,X50,X110,X120,X131,X132,X133;
	int X100,X101,X102,X103,X141,X142,X143,X153,X180,X181,X182,X183,X210,X211;
	int X310, X311, X312, X313, X314, X315, X320, X321, mooring_count, net_count;
	double X21_d,X22_m;
	double X23_x,X23_y,X23_z;
	double X24_Ix,X24_Iy,X24_Iz;	
	double X25_Cp,X25_Cq,X25_Cr;	
    double X26_Cu,X26_Cv,X26_Cw;	
	double X41,X42,X43,X44;
	double X100_x,X100_y,X100_z;
	double X101_phi, X101_theta, X101_psi;
	double X102_u, X102_v, X102_w;
	double X103_p, X103_q, X103_r;
	double *X110_xs,*X110_xe,*X110_ys,*X110_ye,*X110_zs,*X110_ze;
	double X120_rad,X120_xc,X120_yc,X120_zc;
	double X131_rad,X131_h,X131_xc,X131_yc,X131_zc;
	double X132_rad,X132_h,X132_xc,X132_yc,X132_zc;
	double X133_rad,X133_h,X133_xc,X133_yc,X133_zc;
	double X153_xs,X153_xe,X153_ys,X153_ye,X153_zs,X153_ze;
    int X163;
    double *X163_x1,*X163_y1,*X163_z1;
    double *X163_x2,*X163_y2,*X163_z2;
    double *X163_x3,*X163_y3,*X163_z3;
    double *X163_x4,*X163_y4,*X163_z4;
    double *X163_x5,*X163_y5,*X163_z5;
    double *X163_x6,*X163_y6,*X163_z6;
    int X164;
    double *X164_x1,*X164_y1,*X164_z1;
    double *X164_x2,*X164_y2,*X164_z2;
    double *X164_x3,*X164_y3,*X164_z3;
    double *X164_x4,*X164_y4,*X164_z4;
    double *X164_x5,*X164_y5,*X164_z5;
    double *X164_x6,*X164_y6,*X164_z6;
    double *X164_x7,*X164_y7,*X164_z7;
    double *X164_x8,*X164_y8,*X164_z8;
    double X181_x,X181_y,X181_z;
    double X182_x,X182_y,X182_z;
    double X183_x,X183_y,X183_z,X183_phi,X183_theta,X183_psi;
    double X184;
    
    int X205;
    int X206,X207;
    double X206_ts,X206_te,X207_ts,X207_te;
	double X210_u,X210_v,X210_w;
	double X211_p,X211_q,X211_r;
    int X221;
    double X221_xs,X221_xe,X221_ys,X221_ye,X221_zs,X221_ze;
    double *X311_xs,*X311_xe,*X311_ys,*X311_ye,*X311_zs,*X311_ze;
    double *X311_w,*X311_rho_c,*X311_EA,*X311_d,*X311_l,*X311_H,*X311_P,*X311_facT;
    double *X312_k,*X312_T0;
    double *X314_T, *X315_t;
    int *X320_type;
	double *X321_Sn,*X321_d,*X321_lambda,*X321_dk,*X321_rho,*X321_nd,*X321_nl;
    double *X322_D,*X322_L,*X322_x0,*X322_y0,*X322_z0,*X322_phi,*X322_theta,*X322_psi;
    int X324;
    double X323_m,X323_d,X323_l;
    double *X324_x,*X324_y,*X324_z;
    double X325_dt,X325_relX,X325_relY,X325_relZ;
    int X400;
    double X401_p0,X401_cl,X401_cb,X401_a;

    // FSI
    int Z10,Z11,FSI_count;
    double *Z11_x,*Z11_y,*Z11_z,*Z11_l,*Z11_w,*Z11_t,*Z11_rho,*Z11_e,*Z11_ix,*Z11_iy,*Z11_iz,*Z11_nu,*Z11_n;
    double Z12_ckx,Z12_cky,Z12_ckz,Z12_cdx,Z12_cdy,Z12_cdz;
	
	// Grid
	int Y40,Y50,Y60,Y71,Y72,Y73,Y74;

    // Test options
    int Y1,Y2,Y3,Y4,Y5;

	// time + iterations
	int inneriter,count,solveriter,preconiter,count_statestart;
    int solver_status;
    int sediter;
    double final_res;
	double dt,dt_old,simtime,viscmax;
	double mindt,maxdt;
	double umax,vmax,wmax,epsmax,kinmax,pressmin,pressmax,omegamax;
	double presstime,veltime,reinitime,turbtime,plstime,itertime;
	double sedsimtime,sedwavetime;
	double wavetime;
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
    double printtime, sedprinttime,fsfprinttime,probeprinttime,stateprinttime,exportprinttime;
    double partprinttime;

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
	double kintime,epstime;
	double poissontime, laplacetime;
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

	// maxcoor
	double xcoormax,xcoormin,ycoormax,ycoormin,zcoormax,zcoormin;
	double maxlength;


	// wave coefficients
	double wT,wV,wH,wA,wL,wd,ww,wk;
	double wHs,wAs,wwp,ww_s,ww_e,wTp;
	int wN;
    double wts,wte;
    
    // free surface
    double psi,psi0;
	int pressval;

// Boundary
    //int **boundary;
    int **fgc;

	static int knox,knoy,knoz;
	static int margin;

	static int xtp,ytp,ztp;
	static int xmax,ymax,zmax;

// PARALELL
    int mpirank;
	int gcx_1range1[7],gcx_3range1[7];
	int gcx_1range2[7],gcx_3range2[7];
	int gcx_1range3[7],gcx_3range3[7];
	int gcx_1range4[7],gcx_3range4[7];
	
// Non-Uniform Mesh    
    double *XN,*YN,*ZN;
    double *XP,*YP,*ZP;
    double *DXN,*DYN,*DZN;
    double *DXP,*DYP,*DZP;
    double *ZSN,*ZSP;
    double DXM,DYD,DXD;
    double DYM,DZM;
    
    double *RN,*SN,*TN;
    double *RP,*SP,*TP;
    double *DRN,*DSN,*DTN;
    double *DRP,*DSP,*DTP;
    double DRM,DSM,DTM;
    double DX,DY,DZ;
    double *DRDXN,*DSDYN,*DTDZN;
    double *DRDXP,*DSDYP,*DTDZP;
    double *DDRDDXN,*DDSDDYN,*DDTDDZN;
    double *DDRDDXP,*DDSDDYP,*DDTDDZP;
    
    weno_nug_func *wenofunc;
    
// sigma coordinate
    double *sig;
    double *sigx,*sigy,*sigz,*sigt;
    double *sigxx;
    
private:
	void clear(char&, int&);
    
};

#endif
