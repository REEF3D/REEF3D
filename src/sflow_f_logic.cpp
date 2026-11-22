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

#include"sflow_f.h"
#include"sflow_header.h"
#include"sflow_curv.h"
#include"sflow_momentum_SI.h"
#include"sflow_timestep_SI.h"

void sflow_f::logic(lexer *p, fdm2D* b, ghostcell* pgc)
{	
    // curvilinear grid (G 50): options available on it
    if(p->geo_curv==1)
    sflow_curv::check(p);

    // mesh refinement given: a box (G 10), a flag (G 20 surface jump, G 21 second difference, G 22 shoreline,
    // G 23 breaking) or the body zone (G 12).  One test for every choice below that depends on it (eddy viscosity,
    // sflow_amr, the message without anything to refine)
    const bool amr_given = (p->G10>0 || p->G20>0.0 || p->G21>0.0 || p->G22>0 || p->G23==1 || p->G12>0);

    // forcing
    psfdf=new sflow_forcing(p);
    
	// timestep
    if(p->N48==0)
	ptime = new sflow_fixtimestep(p,b);
    
    // N 48 2 (FNPF: stability step of the grid) as N 48 1
    if(p->N48>=1 && p->A226==0)
	ptime = new sflow_etimestep(p,b);
    
    // semi-implicit free surface (A 210 11/12/13): advective time step
    if(p->N48>=1 && p->A226==1)
	ptime = new sflow_timestep_SI(p,b);
	
    // HLL fluxes, signal speeds and reconstruction (A 211: 0 first-order, 1-3 limiter, 4 WENO)
    phll = new sflow_HLL(p,pgc,pBC);
    
    pss = new sflow_signal_speed(p);
    
    if(p->A211<=3)
    precon = new sflow_reconstruct_hires(p,pBC);
    
    if(p->A211>=4)
    precon = new sflow_reconstruct_weno(p,pBC);
    
    // scalar convection (turbulence)
    if(p->A211<=3)
    pconvec = new sflow_fou(p);
    
    if(p->A211>=4)
    pconvec = new sflow_weno_flux(p);

	// free surface
	if(p->A240>=1)
	pfsf = new sflow_eta(p,b,pgc,pBC);

	// diffusion
	if(p->A212==0)
	pdiff =  new sflow_diffusion_void(p);
    
    if(p->A212==0 && p->A260>0 && p->mpirank==0)
    cout<<"SFLOW warning: turbulence model A 260 "<<p->A260<<" with A 212 0 set in ctrl.txt: the eddy viscosity does not enter the momentum equations (leave out A 212 or set 1 or 2)"<<endl;
	
	if(p->A212==1)
	pdiff =  new sflow_ediff(p);
    
    if(p->A212==2)
	pdiff =  new sflow_idiff(p);
	
	// pressure
    if(p->A220==0 || p->A220==4)
	ppress = new sflow_hydrostatic(p,b,pBC);
    
    if(p->A220==1 || p->A220>=5)
	ppress = new sflow_pjm_lin(p,b,pBC);
    
    if(p->A220==2 || p->A220==3)
	ppress = new sflow_pjm_quad(p,b,pgc,pBC);
    
    if(p->A220>=5 && p->mpirank==0)
    cout<<"A 220 "<<p->A220<<" not available with the HLL scheme, using A 220 1"<<endl;

    // A 221 2 (sflow_drywall.h): one Cartesian grid
    if(p->A221==2 && (p->geo_curv==1 || p->geo_mb>0 || p->G1>0) && p->mpirank==0)
    cout<<"SFLOW A 221 2 (dry land next to deep water as a wall): on one Cartesian grid only, here as A 221 1"<<endl;

    // turbulence
	if(p->A260==0)
	pturb =  new sflow_turb_void(p);
    
    if(p->A260==1)
	pturb =  new sflow_turb_ke_IM1(p);
    
    if(p->A260==2)
	pturb =  new sflow_turb_kw_IM1(p);
    
    if(p->A260==3)
	pturb =  new sflow_turb_prandtl(p);
    
    if(p->A260==4)
	pturb =  new sflow_turb_parabolic(p);
    
    // curvilinear grids and REEFAMR (G 1 with boxes, flags or the body zone, multi-block G 50 1): the algebraic eddy
    // viscosity (A 260 3/4) by each grid in its explicit diffusion (sflow_ediff::own_eddyv), from the state of the stage
    const int amr_on = ((p->G1>0 && amr_given) || p->geo_mb>0) ? 1 : 0;
    if((p->geo_curv==1 || amr_on==1) && (p->A260==3 || p->A260==4) && p->A212==1)
    {
    delete pturb;
    pturb = new sflow_turb_void(p);
    static_cast<sflow_ediff*>(pdiff)->own_eddyv = 1;
    }

    if(p->A260==5)
    {
    if(p->mpirank==0)
    cout<<"SFLOW: A 260 5 (k-omega v1) is no longer available, use A 260 2 (depth-averaged k-omega)"<<endl;
    
    exit(1);
    }
    
    // S 16 4: the bed shear stress 0.3 rho k needs the k of A 260 1/2 (zero otherwise: no sediment transport)
    if(p->S10>0 && p->S16==4 && p->A260!=1 && p->A260!=2 && p->mpirank==0)
    cout<<"SFLOW warning: S 16 4 (bed shear stress from k) without the k-epsilon or k-omega model (A 260 1/2): k is zero, no bed shear stress"<<endl;

    // Sediment
    if(p->S10==0)
    psed = new sediment_void();

    if(p->S10>0)
    psed = new sediment_f(p,pgc,nullptr,pBC);     // no CFD turbulence model in SFLOW
	
	// solver
    //  N 10 1: REEFMG, the same multigrid as for FNPF, NHFLOW and CFD, run
    //  on the slice (reefmg2D); otherwise hypre (N 10 11-19)
    if(p->N10==1)
	ppoissonsolv = new reefmg2D(p,pgc);
    
#ifdef REEF3D_USE_HYPRE
    if(p->N10!=1)
	ppoissonsolv = new hypre_struct2D(p,pgc);
#else
    if(p->N10!=1)
	ppoissonsolv = new reefmg2D(p,pgc);
#endif
    
    
    psolv = new sflow_bicgstab(p,pgc);
    
    //IOFlow
	if(p->B60==0 && p->B90==0)
	pflow = new ioflow_v(p,pgc,pBC);

	if(p->B60>=1)
	pflow = new ioflow_f(p,pgc,pBC);

	if(p->B90>=1)
	pflow= new iowave(p,pgc,pBC);

    // B 60 3/4: the outflow hydrograph is read and printed (Qo_ipol) but not applied: the outflow of SFLOW is the
    // water level F 62 (F 50 2/3) or free; B 60 3 has no inflow either (W 10 with B 60 1, hydrograph.dat with B 60 2/4)
    if(p->B60>=3 && p->B90==0 && p->mpirank==0)
    cout<<"SFLOW warning: B 60 "<<p->B60<<": the outflow hydrograph is not applied in SFLOW (the outflow is the water level F 62 with F 50 2)"
        <<(p->B60==3 ? ", and B 60 3 has no inflow (B 60 1: W 10, B 60 2: hydrograph.dat)" : "")<<endl;

    // B 60 with the explicit integrators and no water level at the outflow: the ghost cells are zero gradient in the
    // level and the velocity, nothing sets the level of a subcritical river (a channel with the normal depth
    // 1.47 m filled to 4.8 m); the semi-implicit step has a free outflow there (Riemann invariant)
    if(p->B60>=1 && p->B90==0 && p->A226==0 && p->F50!=2 && p->F50!=3 && p->B440==0 && p->B441==0 && p->B442==0 && p->B520==0 && p->mpirank==0)
    cout<<"SFLOW warning: B 60 "<<p->B60<<" with A 210 "<<p->A210<<" (explicit) and F 50 "<<p->F50<<": no water level at the outflow, the level of a subcritical flow is not determined (set F 50 2 and F 62, or A 210 11 for a free outflow)"<<endl;

	
	// printer
	pprint = new sflow_vtp_fsf(p,b,pgc);
	
	pprintbed = new sflow_vtp_bed(p);
    
    //6DOF
    if(p->X10<2)
    p6dof = new sixdof_void(p,pgc);
    
    if(p->X10>=2)
    p6dof = new sixdof_sflow(p,pgc);
	
	// momentum
    if(p->A210==2 && p->A226==0)
	pmom = new sflow_momentum_RK2(p,b,pgc,phll,pss,precon,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,psfdf,p6dof);
    
	if(p->A210!=2 && p->A226==0)
	pmom = new sflow_momentum_RK3(p,b,pgc,phll,pss,precon,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,psfdf,p6dof);
    
    // semi-implicit free surface, theta method (A 210 11/12/13, A 227), separate from the HLL integrators
    if(p->A226==1)
	pmom = new sflow_momentum_SI(p,b,pgc,phll,pss,precon,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,psfdf,p6dof);

    // mesh refinement (G 1): static boxes (G 10) and/or flagged regions (G 20-23)
    if(amr_on==1)
    {
    sflow_momentum_RK3 *prk3 = dynamic_cast<sflow_momentum_RK3*>(pmom);
    sflow_momentum_SI *psi = dynamic_cast<sflow_momentum_SI*>(pmom);

        // semi-implicit free surface (A 210 12/13): the blocks and the patches have their own semi-implicit objects; an RK3
        // object of level 0 serves the fills and the restriction of the refinement (sflow_amr_si.cpp)
        if(psi!=nullptr)
        {
        met_forcing *pm = b->met;       // the met forcing of the semi-implicit object (the helper makes its own)
        sflow_momentum_RK3 *hrk3 = new sflow_momentum_RK3(p,b,pgc,phll,pss,precon,pdiff,ppress,psolv,ppoissonsolv,pflow,pfsf,psfdf,p6dof);
        b->met = pm;
        pamr = new sflow_amr(p,b,pgc,pBC,p6dof,pflow,phll,hrk3);
        pamr->si_attach(psi);
        }

        if(prk3!=nullptr)
        {
        pamr = new sflow_amr(p,b,pgc,pBC,p6dof,pflow,phll,prk3);
        prk3->pamr = pamr;
        
        // non-hydrostatic pressure on the composite grid
        sflow_pressure_nh *pnh = dynamic_cast<sflow_pressure_nh*>(ppress);
        if(p->A220>=1 && p->A220<=3 && pnh!=nullptr)
        {
        pnh->amr = pamr;
        pamr->pnh0 = pnh;
        }
        }

        // the refinement runs with the RK3 integrator only: say so instead of running the native grid
        // alone (RK2 A 210 2, the semi-implicit free surface A 210 11/12/13)
        if(prk3==nullptr && psi==nullptr && p->G1>0 && p->geo_mb==0 && p->mpirank==0)
        cout<<"SFLOW AMR (G 1): only with A 210 3, 11, 12 or 13 -- refinement switched off"<<endl;

        // multi-block grid: the blocks 1.. run in the mesh refinement (RK3 or the semi-implicit free surface)
        if(prk3==nullptr && psi==nullptr && p->geo_mb>0)
        {
            if(p->mpirank==0)
            cout<<"!!! SFLOW multi-block grid (G 50 1 with DIVEMesh R 20 1): only with A 210 3, 11, 12 or 13 !!!"<<endl;
            exit(1);
        }

        // A 260 1/2: the level-0 turbulence model for the models of the patches (sflow_amr_turb.cpp)
        if(pamr!=nullptr)
        pamr->turb_attach(pturb);
    }

    // G 1 without anything to refine: no box (G 10), no flag (G 20-23) and no body zone (G 12)
    if(p->G1>0 && !amr_given && p->geo_mb==0 && p->mpirank==0)
    cout<<"SFLOW AMR (G 1): no refinement box (G 10), flag (G 20-23) or body zone (G 12) given -- refinement switched off"<<endl;
    
    // two-layer landslide (W 150, slide2D): the explicit HLL (A 210 2/3) and the semi-implicit (A 210 11/12/13)
    // water on one Cartesian grid; the semi-implicit step keeps the water column of the end of the step before
    // (si_restore) and takes the moved bed through depth_update, as with sediment
    if(p->W150>0)
    {
        string why;

        if(p->S10>0)
        why = "sediment transport (S 10): both move the bed";

        if(p->geo_curv==1)
        why = "curvilinear grids (G 50)";

        if(pamr!=nullptr || p->geo_mb>0)
        why = "mesh refinement and multi-block grids (G 1, G 50 1)";

        // the quadratic non-hydrostatic pressure (A 220 2) takes the moving bed in its bed condition; the others not
        if(p->A220>0 && p->A220!=2)
        why = "the non-hydrostatic pressure or Boussinesq A 220 "+to_string(p->A220)+": A 220 0 or 2";
        
        if(p->A220==2 && p->A226==1)
        why = "the non-hydrostatic pressure A 220 2 and the semi-implicit water A 210 "+to_string(10+p->A210)+" together (not tested)";

        if(!why.empty())
        {
            if(p->mpirank==0)
            cout<<"!!! SFLOW slide layer (W 150): not available with "<<why<<" !!!"<<endl;
            pgc->final(true);
        }

    }

    // F 74: the staged release of the material of F 73 (sflow_release, sflow_level_file.cpp)
    if(p->F74>0)
    {
        string why;

        if(p->F73==0)
        why = "without F 73 (the material comes from sflow-level.dat)";

        if(p->W150>0)
        why = "the slide layer (W 150)";

        if(p->S10>0)
        why = "sediment transport (S 10): both move the bed";

        if(pamr!=nullptr || p->geo_mb>0)
        why = "mesh refinement and multi-block grids (G 1, G 50 1)";

        if(p->I45==1)
        why = "a restart from a checkpoint (I 45 1): the release state is not in the checkpoints yet";

        if(!why.empty())
        {
            if(p->mpirank==0)
            cout<<"!!! SFLOW staged release (F 74): not available with "<<why<<" !!!"<<endl;
            pgc->final(true);
        }
    }

    //Potential Flow Solver
    if(p->I11==0)
    potflow = new sflow_potential_v(p);

    if(p->I11==1)
    potflow = new sflow_potential_f(p);
    
}
