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

#include"nhflow_curv.h"
#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"vrans_header.h"
#include"nhflow_header.h"
#include"6DOF_void.h"
#include"6DOF_nhflow.h"
#include"dem_header.h"
#include"nhflow_amr.h"
#include"nhflow_longrun.h"
#include<fstream>
#include<sstream>
#include<string>

// membranes (X 330): true if membrane.dat contains a moving membrane (structure rigid or flexible); rank 0 reads
static bool nhflow_membrane_moving(lexer *p, ghostcell *pgc)
{
    int moving=0;

    if(p->mpirank==0)
    {
        std::ifstream f("membrane.dat");
        std::string line;

        while(std::getline(f,line))
        {
            std::istringstream ls(line.substr(0,line.find('#')));
            std::string key, val;

            if(ls>>key>>val && key=="structure" && (val=="rigid" || val=="flexible"))
            moving=1;
        }
    }

    MPI_Bcast(&moving,1,MPI_INT,0,pgc->mpi_comm);

    return moving==1;
}

void driver::logic_nhflow()
{    
	if(p->mpirank==0)
    cout<<"creating objects"<<endl;
    
    p->phimean = p->wd = p->F60;
    
    // two-layer landslide (W 150, slide2D): the explicit and the semi-implicit NHFLOW on one Cartesian grid
    if(p->W150>0)
    {
        std::string why;
        
        if(p->S10>0)
        why = "sediment transport (S 10): both move the bed";
        
        if(p->geo_curv==1)
        why = "curvilinear grids (G 50)";
        
        if(p->G1>0 || p->geo_mb>0)
        why = "mesh refinement and multi-block grids (G 1, G 50 1)";
        
        if(!why.empty())
        {
            if(p->mpirank==0)
            cout<<"!!! NHFLOW slide layer (W 150): not available with "<<why<<" !!!"<<endl;
            pgc->final(true);
        }
    }
    
    // river free-surface initialisation (F 65): reference level at the outflow when F 60 is not given
    if(p->F65>0 && p->F60<-1.0e19 && p->F62>-1.0e19)
    p->phimean = p->wd = p->F62;
    
    // curvilinear grid (G 50): unsupported options are refused
    if(p->geo_curv==1)
    nhflow_curv::check(p);
    
// nhflow
    if(p->A10!=5)
    pnhf=new nhflow_v(p,d,pgc);
    
    if(p->A10==5)
    pnhf=new nhflow_f(p,d,pgc);
    
// forcing
    pnhfdf=new nhflow_forcing(p,d,pgc);

// FSF
    pnhfsf = new nhflow_fsf_f(p,d,pgc,pflow,pBC);
    
// time stepping
    // time stepping; hydrostatic semi-implicit free surface (A 510 11/12/13): the advective time step
    if(p->A534==1)
    {
    nhflow_momentum_SI::check(p);
    pnhfstep=new nhflow_timestep_SI(p);
    }
    
    if(p->A534!=1)
	pnhfstep=new nhflow_timestep(p);

    // density from the tracers (A 620): needs the tracers
    if(p->A620>0 && p->A600==0)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW density (A 620) needs the tracers (A 600 1/2) !!!"<<endl;
        exit(1);
    }

    // mixing damped by the stratification (A 623): needs the density; A 547 or k-epsilon
    if(p->A623>0 && p->A620==0)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW stratified mixing (A 623) needs the density (A 620 1/2) !!!"<<endl;
        exit(1);
    }

    if(p->A623>0 && p->A560!=0 && p->A560!=1 && p->A560!=21)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW stratified mixing (A 623): with the algebraic eddy viscosity (A 547) or k-epsilon (A 560 1, 21) only !!!"<<endl;
        exit(1);
    }

    // semi-Lagrangian momentum advection (A 549 1): the hydrostatic semi-implicit path
    if(p->A549!=0 && (p->A549!=1 || p->A534!=1 || p->A549_c<=0.0 || p->A549_c>2.5))
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW A 549 m C: m 0 or 1 (semi-Lagrangian momentum advection, with A 510 11/12/13 only), 0 < C <= 2.5 !!!"<<endl;
        exit(1);
    }

    // eddy viscosity floor of a turbulence model in thin water (A 627): the semi-implicit path with a turbulence model
    if(p->A627_c>0.0 && (p->A534!=1 || p->A560==0))
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW eddy viscosity floor (A 627): with the hydrostatic semi-implicit free surface (A 510 11/12/13) and a turbulence model (A 560 > 0) only !!!"<<endl;
        exit(1);
    }

    // tracers (A 600): the semi-implicit path only
    if(p->A600>0 && p->A534!=1)
    {
        if(p->mpirank==0)
        cout<<"!!! NHFLOW tracers (A 600) need the hydrostatic semi-implicit free surface (A 510 11/12/13) !!!"<<endl;
        exit(1);
    }

//discretization scheme
    // signal speed
    pss = new nhflow_signal_speed(p);
    
    // A 520: 0 hydrostatic, 1 projection, 2 incremental projection (A 520 3, the HLLYL flux with its pressure
    // scheme nhflow_pjm_yl, has been removed)
    if(p->A520<0 || p->A520>2)
    {
        if(p->mpirank==0)
        cout<<endl<<"A 520 "<<p->A520<<": unknown NHFLOW pressure scheme; use 0 (hydrostatic), 1 (projection) or 2 (incremental projection). "
            <<"A 520 3 (the HLLYL flux) has been removed."<<endl<<endl;
        
        pgc->final(true);
    }
    
//Convection	
    if(p->A511==1)
	pnhfconvec = new nhflow_HLL(p,pgc,pBC);
    
    if(p->A511==2)
	pnhfconvec = new nhflow_HLL(p,pgc,pBC);   // HLLC mode inside nhflow_HLL (A 511 2)
    
    pnhfscalarconvec = new nhflow_scalar_ifou(p,1);   // advective form for k, epsilon/omega
    
//Diffusion
    if(p->A512==0 && p->A560==0)
    pnhfdiff = new nhflow_diff_void(p);
    
    if(p->A512==0 && p->A560>0 && p->j_dir==1)
    pnhfdiff = new nhflow_idiff(p);
    
    if(p->A512==0 && p->A560>0 && p->j_dir==0)
    pnhfdiff = new nhflow_idiff_2D(p);
    
    if(p->A512==1)
    pnhfdiff = new nhflow_ediff(p);
    
    if(p->A512==2 && p->j_dir==1)
    pnhfdiff = new nhflow_idiff(p);
    
    if(p->A512==2 && p->j_dir==0)
    pnhfdiff = new nhflow_idiff_2D(p);
    
// reconstruction    
    if(p->A514<=3)
    precon = new nhflow_reconstruct_hires(p,pBC);
    
    if(p->A514==4 || p->A514==5)
    precon = new nhflow_reconstruct_weno(p,pBC);
    
//pressure scheme
    // moving membranes (X 330, structure rigid or flexible) need the full projection: with the incremental scheme
    // the old pressure in the membrane layer, where the mobility is ~ 1/(1 + a K) << 1, is only weakly controlled by
    // the projection and drifts; when the layer moves, cells leaving it release that pressure (2D flexible bag:
    // unstable after 2.7 s, collar cases within 2 s)
    if(p->X330>0 && p->A520==2 && nhflow_membrane_moving(p,pgc))
    {
        if(p->mpirank==0)
        cout<<"X 330: moving membrane, pressure scheme A 520 2 -> A 520 1 (full projection)"<<endl;

        p->A520=1;
    }
    
    // the hull as a sharp thin body (X 15 8) keeps the incremental projection (nhflow_pjm_corr: the old pressure in the
    // predictor with the link mobilities, the wall velocity in the right-hand side). The full projection damps short
    // waves: a 2D wave of T 1.4 s in 3 m depth (dx 0.1 m, 30 layers) kept 88.5 % of its height over 20 m with A 520 1
    // against 98.5 % with A 520 2; in the ship frame of the KVLCC2-like hull (dx 0.3 m, 20 layers) the incident wave
    // reached the hull at 62 % of its height at lambda/L 0.6 and 93 % at 1.0 (Ship Hydrodynamics validation 34).
    // On the mesh refinement hierarchy the hull needs the full projection (nhflow_amr::keys_ok).
    if(p->X15==8 && p->X10>0 && p->A520==2 && p->G1>0)
    {
        if(p->mpirank==0)
        cout<<"X 15 8 on the mesh refinement hierarchy (G 1): pressure scheme A 520 2 -> A 520 1 (full projection)"<<endl;

        p->A520=1;
    }
    
    // the body as a solid with cut faces (X 15 10) keeps the incremental projection too (nhflow_pjm_corr: the
    // finite-volume pressure force of the cut cells in the predictor and the correction, the pressure extension into
    // the nodes without water; Ship Hydrodynamics validation 43). On the mesh refinement hierarchy it needs the full
    // projection (nhflow_amr::keys_ok).
    if(p->X15==10 && p->X10>0 && p->A520==2 && p->G1>0)
    {
        if(p->mpirank==0)
        cout<<"X 15 10 on the mesh refinement hierarchy (G 1): pressure scheme A 520 2 -> A 520 1 (full projection)"<<endl;

        p->A520=1;
    }

    // direct forcing with the body in the projection (X 15 11): needs the projection, on the uniform grid; else the
    // forcing of X 15 4
    if(p->X15==11 && p->X10>0 && (p->A520==0 || p->G1>0 || p->geo_curv==1))
    {
        if(p->mpirank==0)
        cout<<"X 15 11 is not available "<<(p->A520==0 ? "without the pressure projection (A 520 0)" : (p->G1>0 ? "on the mesh refinement hierarchy (G 1)" : "on the curvilinear grid (G 50)"))<<" yet: X 15 4"<<endl;

        p->X15=4;
    }

    if(p->A520==0)
	pnhpress = new nhflow_pjm_hs(p,d,pBC);
    
    if(p->A520==1)
    pnhpress = new nhflow_pjm(p,d,pgc,pBC);
    
    if(p->A520==2)
    pnhpress = new nhflow_pjm_corr(p,d,pgc,pBC);
    
    // A 508 1 with the explicit schemes: fourth-order horizontal operators of the projection (nhflow_pjm_ho.h)
    if(p->A508==1 && p->A534==0 && p->A520>0 && (p->geo_curv==1 || p->G1>0 || p->geo_mb>0 || p->X10>0) && p->mpirank==0)
    cout<<"NHFLOW: A 508 1 (fourth-order horizontal operators of the projection) needs one uniform Cartesian grid without 6DOF, second order used"<<endl;

//Turbulence
    // unknown A 560: no model would be created (pnhfturb uninitialised) and the run crashed later
    if(p->A560!=0 && p->A560!=1 && p->A560!=21 && p->A560!=2 && p->A560!=22 && p->A560!=31)
    {
        if(p->mpirank==0)
        cout<<endl<<"A 560: unknown NHFLOW turbulence model; use 0 (laminar), 1/21 (k-epsilon), 2/22 (k-omega) or 31 (LES Smagorinsky)."<<endl<<endl;
        
        pgc->final(true);
    }
    
    // no bed friction in the momentum equations (A 519 0; the diffusion has a zero-gradient bed with A 513 1): the flow
    // slips over the bed, a RANS model gets its turbulence only from the wall function of the bottom cell and the
    // velocity profile stays a plug (open-channel flow: nu_t far below the parabolic profile, the bed shear stress of
    // the sediment S 16 from the near-bed velocity about twice too high)
    if(p->A519==0 && p->A518!=2 && (p->A560==1 || p->A560==21 || p->A560==2 || p->A560==22 || p->S10>0) && p->mpirank==0)
    cout<<"NOTE A 519 0: no bed friction in the momentum equations, the flow slips over the bed; A 519 1 (the default "
        <<"with a turbulence model A 560 or sediment transport S 10) for open-channel flow"<<endl;
    
    // k-epsilon and k-omega without the stabilisation of Larsen & Fuhrman (2018, A 565 1) in waves: the strain of the
    // (nearly) potential wave motion produces k in the whole water column with nothing to dissipate it, and nu_t
    // grows without bound. Regular wave H 0.06 m, T 1.5 s in 0.5 m of water, 30 m flume, 10 layers: after 200
    // periods nu_t = 5e-3 m2/s in the whole column and the wave height at 23 m down to 45 % (k-epsilon and k-omega,
    // A 564 1); with A 565 1 nu_t stays at 3e-6 m2/s outside the bed layer and the wave height at 99 %
    if(p->B90>0 && p->A565==0 && (p->A560==1 || p->A560==21 || p->A560==2 || p->A560==22) && p->mpirank==0)
    cout<<"NOTE A 565 0 with waves (B 90): the eddy viscosity of k-epsilon / k-omega grows without bound under "
        <<"non-breaking waves and damps them; A 565 1 (stabilisation of Larsen & Fuhrman 2018) for wave runs"<<endl;
    
    // bed of the velocity (A 518): 1 the wall-function bed, mirrored ghost cells, zero gradient in the diffusion and
    // the bed shear from the log-law friction A 519; 2 the resolved no-slip bed, u = 0 on the bed face in the
    // ghost cells and in the diffusion, the bed shear from the diffusion (gc_startV.cpp, nhflow_idiff, nhflow_bcmom)
    if(p->A518==2 && p->mpirank==0)
    {
        if(p->A519>0)
        cout<<"NOTE A 518 2: resolved no-slip bed, the bed shear comes from the diffusion; no bed friction of A 519 at the bed"<<endl;
        
        if(p->A512==0 && p->A560==0)
        cout<<"!!! A 518 2 without diffusion (A 512 0, A 560 0): the no-slip bed acts only through the ghost cells"<<endl;
        
        // the bed face has the molecular viscosity only (nu_t = 0 at the wall): with a RANS model the bed cell has to lie
        // in the viscous sublayer (z+ < 5); on a coarse grid the bed shear is far too low (open channel, k-epsilon,
        // 20 layers over 2 m: u* half the log law), A 518 1 there
        if(p->A560==1 || p->A560==21 || p->A560==2 || p->A560==22)
        cout<<"!!! A 518 2 with a RANS model (A 560 "<<p->A560<<"): wall-resolving grids only (bed cell z+ < 5), A 518 1 otherwise"<<endl;
        
        if(p->A534==1)
        cout<<"!!! A 518 2 with the semi-implicit free surface (A 510 11/12/13): its column operator keeps the log-law bed friction (A 519 1)"<<endl;
    }
    
    // HLES split (A 561) and QCR (A 563) work with the RANS models (k-epsilon, k-omega) only
    if((p->A561!=0 || p->A563>0.0) && p->A560!=1 && p->A560!=21 && p->A560!=2 && p->A560!=22)
    {
        if(p->mpirank==0)
        cout<<"A 561 / A 563: need a RANS model (A 560 1, 21, 2 or 22), switched off"<<endl;
        
        p->A561=0;
        p->A563=0.0;
    }
    
    // EARSM (A 559): k-omega only; it replaces the QCR (A 563)
    if(p->A559!=0 && p->A560!=2 && p->A560!=22)
    {
        if(p->mpirank==0)
        cout<<"A 559: EARSM needs k-omega (A 560 2 or 22), switched off"<<endl;
        
        p->A559=0;
    }
    
    if(p->A559!=0 && p->A563>0.0)
    {
        if(p->mpirank==0)
        cout<<"A 559: EARSM on, the QCR (A 563) is switched off"<<endl;
        
        p->A563=0.0;
    }
    
    if(p->A561<0 || p->A561>2)
    {
        if(p->mpirank==0)
        cout<<endl<<"A 561: unknown HLES option; use 0 (off), 1 (nu_h = nu_t + nu_sgs) or 2 (nu_h = nu_sgs)."<<endl<<endl;
        
        pgc->final(true);
    }
    
    if(p->A560==0)
	pnhfturb = new nhflow_komega_func_void(p,d,pgc);
    
    if(p->A560==1 || p->A560==21)
    {
	pnhfturb = new nhflow_kepsilon_IM1(p,d,pgc);
    
        if(p->j_dir==1)
        pnhfturbdiff = new nhflow_idiff(p);
        
        if(p->j_dir==0)
        pnhfturbdiff = new nhflow_idiff_2D(p);
    }
    
    if(p->A560==2 || p->A560==22)
    {
	pnhfturb = new nhflow_komega_IM1(p,d,pgc);
    
        if(p->j_dir==1)
        pnhfturbdiff = new nhflow_idiff(p);
        
        if(p->j_dir==0)
        pnhfturbdiff = new nhflow_idiff_2D(p);
    }
    
    if(p->A560==31)
	pnhfturb = new nhflow_LES_Smagorinsky(p,d,pgc);

//Solver
    if(p->j_dir==0)
	psolv = new bicgstab_ijk_2D(p,a,pgc);
    
    if(p->j_dir==1)
	psolv = new bicgstab_ijk(p,a,pgc);

//Poison Solver	
	if(p->N10==0)
	ppoissonsolv = new solver_void(p,a,pgc);
    
    if(p->N10==1)
    ppoissonsolv = new reefmg(p,pgc,p->N10,p->N11);
    
    if(p->N10==3 && p->j_dir==0)
	ppoissonsolv = new bicgstab_ijk_2D(p,a,pgc);
    
    if(p->N10==3 && p->j_dir==1)
	ppoissonsolv = new bicgstab_ijk(p,a,pgc);
	
#ifdef REEF3D_USE_HYPRE
	if(p->N10>=10 && p->N10<20)
	ppoissonsolv = new hypre_struct(p,pgc,p->N10,p->N11);
    
	if(p->N10>=20 && p->N10<30)
	ppoissonsolv = new hypre_aij(p,a,pgc);
    
	if(p->N10>=30 && p->N10<40)
	ppoissonsolv = new hypre_sstruct(p,a,pgc);
#endif

//Printer
    if(p->P150==0)
	pdata = new expdata_void(p,a,pgc);

	if(p->P150>0)
	pdata = new expdata_f(p,a,pgc);
    
    pprint = new printer_nhflow(p,d,pgc);

//VRANS
    if(p->B200==0)
	pnhfvrans = new vrans_nhflow_v(p,d,pgc);

	if(p->B200==1)
	pnhfvrans = new vrans_nhflow_f(p,d,pgc);

//IOFlow
    // the zones of the native grid (B 520 without B 522), their backgrounds (B 510-515), the wave sources (B 500-504)
    // and the start with the background current (B 531) belong to iowave (B 90 1): without it they are read and not
    // used. The zones on the faces of blocks 1.. (B 522) work without it (nhflow_amr, ioflow_amr_block)
    if(p->B90==0 && p->mpirank==0)
    {
        const int native = p->B520 - p->B522;
        const int bg = p->B510 + p->B511 + p->B513 + p->B514 + p->B515;
        const int src = p->B500 + p->B501 + p->B502 + p->B504;
        std::string what;
        auto add = [&](const std::string &w) { what += (what.empty() ? "" : ", ") + w; };
        if(native>0)
        add(std::to_string(native)+" zone(s) of B 520 on the native grid");
        if(bg>0 && p->B522==0)
        add("the backgrounds of B 510-515");
        if(src>0)
        add("the wave sources of B 500-504");
        if(p->B531!=0)
        add("B 531 (the start with the background current)");
        if(!what.empty())
        cout<<"!!! B 90 0: "<<what<<" given but not used: they need the wave module, B 90 1 (a current without waves: B 90 1 with waves of zero height) !!!"<<endl;
    }

	if(p->B60==0 && p->B90==0 && p->B180==0)
	pflow = new ioflow_v(p,pgc,pBC);

	if(p->B60>=1)
	pflow = new ioflow_f(p,pgc,pBC);

	if(p->B90>=1)
	pflow = new iowave(p,pgc,pBC);

	if(p->B180==1||p->B191==1||p->B192==1)
	pflow = new ioflow_gravity(p,pgc,pBC);
    
//Potential Flow Solver
    if(p->I11==0)
    pnhfpot = new nhflow_potential_v;

    if(p->I11==1)
    pnhfpot = new nhflow_potential_f(p);
    
//6DOF
    if(p->X10==0)
    p6dof = new sixdof_void(p,pgc);
    
    if(p->X10>0)
    p6dof = new sixdof_nhflow(p,pgc);

// DEM
    if(p->E10>0)
    pdem = new dem_f(p,pgc);
    else
    pdem = new dem_void();

    pnhfdf->pdem = pdem;
    
// Sediment
    // S 10 for NHFLOW: 1 the sediment step every S 44 flow steps (every momentum solver), 12 the bed advanced in both
    // stages of the SSP-RK2 momentum (A 510 2 only: no other momentum solver calls its stages, the bed never moved and
    // nothing was said); any other S 10 left the sediment object uninitialised. S 12 0 or 1; the bed shear stress
    // formulations of bedshear_nhflow (any other S 16 returned zero shear stress without a message)
    {
    int stop=0;
    
        if(p->S10!=0 && p->S10!=1 && p->S10!=12)
        {
        if(p->mpirank==0)
        cout<<endl<<"!!! S 10 "<<p->S10<<": NHFLOW has S 10 1 (sediment step every S 44 flow steps) and S 10 12 (bed in both stages of the RK2 momentum, A 510 2) !!!"<<endl<<endl;
        stop=1;
        }
        
        if(p->S10==12 && (p->A510!=2 || p->A534==1))
        {
        if(p->mpirank==0)
        cout<<endl<<"!!! S 10 12 (bed advanced in both stages of the RK2 momentum) needs A 510 2; use S 10 1 with "
            <<(p->A534==1 ? "the semi-implicit momentum (A 510 11, 12, 13)" : "A 510 3")<<" !!!"<<endl<<endl;
        stop=1;
        }
        
        if(p->S10>0 && p->S12!=0 && p->S12!=1)
        {
        if(p->mpirank==0)
        cout<<endl<<"!!! S 12 "<<p->S12<<": suspended load is S 12 0 (off) or 1 (van Rijn reference concentration) !!!"<<endl<<endl;
        stop=1;
        }
        
        if(p->S10>0 && !(p->S16==1 || p->S16==2 || p->S16==3 || p->S16==4 || p->S16==6 || p->S16==9 || p->S16==10 || p->S16==11))
        {
        if(p->mpirank==0)
        cout<<endl<<"!!! S 16 "<<p->S16<<": the bed shear stress of NHFLOW is S 16 1, 2, 3, 4, 6, 9, 10 or 11 !!!"<<endl<<endl;
        stop=1;
        }
        
        if(stop==1)
        pgc->final(true);
    
    // S 50 (fixed bed at the inflow / outflow) acts on the CFD bed (sediment_exner, the reinitialisation of the topography)
    if(p->S10>0 && p->S50!=4 && p->mpirank==0)
    cout<<"NOTE S 50 "<<p->S50<<": the bed boundary condition S 50 is CFD only; the edges of the NHFLOW bed are set by the "
        <<"inflow guard S 78 and the sediment inflow S 129 / S 130 (fixed bed zones: S 76, no erosion: S 74, S 134 1 for solids)"<<endl;
    
    // S 12 without a turbulence model: the diffusivity of the suspended load is the molecular viscosity
    if(p->S10>0 && p->S12>0 && p->A560==0 && p->mpirank==0)
    cout<<"NOTE S 12 with A 560 0: without a turbulence model the suspended load has no eddy diffusivity (molecular "
        <<"viscosity only), the sediment stays in the bottom layer; A 560 1, 2 for the mixing over the depth"<<endl;
    }
    
    if(p->S10==0)
    psed = new sediment_void();

    if(p->S10==1)
    psed = new sediment_f(p,pgc,pturbcfd,pBC);
    
    if(p->S10==12)
    psed = new sediment_RK2(p,pgc,pturbcfd,pBC);
      
//Momentum
    if(p->A510==2 && p->A534!=1)
	pnhfmom = new nhflow_momentum_RK2(p,d,pgc,p6dof,pnhfvrans,pnhfdf,psed);

    if(p->A510==3 && p->A534!=1)
	pnhfmom = new nhflow_momentum_RK3(p,d,pgc,p6dof,pnhfvrans,pnhfdf);    
    
    // hydrostatic semi-implicit free surface (A 510 11/12/13, theta method)
    if(p->A534==1)
    {
    nhflow_momentum_SI *psi = new nhflow_momentum_SI(p,d,pgc,p6dof,pnhfvrans,pnhfdf);
    psi->pBC = pBC;
	pnhfmom = psi;
    }
    

// mesh refinement (G 1): patches on the NHFLOW grid
    // multi-block grid (G 50 1, DIVEMesh R 20 1): the blocks 1.. are level-0 patches of the hierarchy
    if(p->G1>0 || p->geo_mb>0)
    pnhfamr = new nhflow_amr(p,d,pgc,pnhfmom,pnhfconvec,pnhfstep,p6dof);

// long-run statistics (P 160): harmonic analysis and time averages
    if(p->P160>0)
    plongrun = new nhflow_longrun(p);

//Lagrangian particles
    if(p->L10==0)
    pnhfpart = new nhflow_particle_void();
    if(p->L10>0)
    pnhfpart = new nhflow_particle_f(p,pgc);
}