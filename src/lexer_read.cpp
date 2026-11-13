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

#include"lexer.h"
#include"ghostcell.h"

void lexer::lexer_read(ghostcell *pgc)
{
    ii_recv=ii_send=0;
    dd_recv=dd_send=0;

    if(mpirank==0)
    control::read_control(this);

    // SFLOW: A 212 not given in ctrl.txt -> explicit horizontal diffusion when a turbulence model is on (A 260 > 0),
    // so that the eddy viscosity reaches the momentum equations; no diffusion otherwise (the previous default)
    if(mpirank==0 && A212<0)
    A212 = (A260>0) ? 1 : ((A246==4) ? 2 : 0);   // A 246 4: the breaking viscosity acts through the implicit diffusion (explicit: the viscosity grows within the step beyond the limit of the step)

    // FNPF: A 316 not given -> the first-order wet-dry gradient fallback only with the moving wet-dry front
    // (A 343 2/3); with the static coastline (A 343 0/1) the legacy zero gradient: the fallback lowered the
    // wave heights at the gauges of static-coastline cases and made them unstable without coastline viscosity
    if(mpirank==0 && A316<0)
    A316 = (A343>=2) ? 1 : 0;

    // SFLOW: A 210 11 / 12 / 13 the semi-implicit theta step (sflow_momentum_SI, A 227, A 228, wetting and drying A 243) with a
    // one-step (BCG advection) / SSP-RK2 / RK3 predictor of the explicit terms; inside REEF3D A226 = 1 marks the semi-implicit
    // step and A210 = 1 / 2 / 3 is its predictor (A 226 itself is no longer a ctrl.txt key)
    if(mpirank==0 && (A210==11 || A210==12 || A210==13))
    {
    A226 = 1;
    A210 -= 10;
    }

    // SFLOW: any other A 210 would run the explicit RK3 step without a word (sflow_f_logic: A210 != 2 -> RK3)
    if(mpirank==0 && A10==2 && A226==0 && A210!=2 && A210!=3)
    {
    cout<<"!!! ctrl.txt: A 210 "<<A210<<" is not a SFLOW time scheme: explicit A 210 2 / 3 (SSP-RK2 / RK3), semi-implicit A 210 11 / 12 / 13 (one-step / SSP-RK2 / SSP-RK3 predictor) !!!"<<endl;
    exit(1);
    }

    // N 47 not given in ctrl.txt -> 0.3; SFLOW A 210 11 (inside A226 = 1, A210 = 1): 0.5, the flow Courant number at most
    // 0.4 there (sflow_timestep_SI: u dt/dx <= 0.8 N 47), the accuracy limit of the one-step predictor at bores; with 0.3 a
    // flow-limited run took 1.5 - 2.6 times as long as A 210 13 for the same time error
    if(mpirank==0 && N47<0.0)
    N47 = (A10==2 && A226==1 && A210==1) ? 0.5 : 0.3;

    // SFLOW: wetting and drying is A 243 for both schemes; the semi-implicit step (A 210 12/13) runs its Newton
    // iteration on the wet pattern (Casulli 2009) for A 243 > 0 (inside REEF3D: A229 = 1; A 229 itself is no longer a
    // ctrl.txt key)
    if(mpirank==0)
    A229 = (A243>0) ? 1 : 0;

    // NHFLOW: A 510 11 / 12 / 13 the semi-implicit theta step (nhflow_momentum_SI, A 535, A 536) with a one-step (BCG
    // advection) / SSP-RK2 / SSP-RK3 predictor of the explicit terms; inside REEF3D A534 = 1 marks the
    // semi-implicit step and A510 = 1 / 2 / 3 its predictor (A 534 itself is no longer a ctrl.txt key)
    if(mpirank==0 && (A510==11 || A510==12 || A510==13))
    {
    A534 = 1;
    A510 -= 10;
    }

    // CFD: W 130 > 0 the rheology core; inside REEF3D W90 > 0 selects the rheology fluid update, density and
    // pressure classes of the CFD logic (W 90 itself stays the legacy model when W 130 is 0)
    if(mpirank==0 && A10==6 && W130>0)
    {
        if(W90>0)
        cout<<"W 130 "<<W130<<" replaces W 90 "<<W90<<endl;
        
    W90 = 1;
    W110 = 1;
    
    // the implicit momentum diffusion with viscosity contrasts of 1e4 - 1e7 needs many more solver iterations than
    // the default cap: with N 46 250 the u, v, w systems stop unconverged and the yield surfaces oscillate
        if(N46<2000)
        {
        cout<<"W 130: N 46 "<<N46<<" -> 2000 (solver iterations of the implicit momentum diffusion with the rheology)"<<endl;
        N46 = 2000;
        }

    // the wall boundary of the implicit diffusion: D 22 1 folds the wall ghost cells into the diagonal (zero
    // gradient), which leaves the wall stress to the log law of bcmom; a viscous or viscoplastic material has no wall
    // law and slipped freely along the walls (viscous gravity current, case 09 of Validation/rheology: front 0.31 m at
    // 1 s against 0.14 m of Huppert 1982). D 22 2 takes the ghost cells of the boundary conditions (no slip)
        if(D22==1)
        {
        cout<<"W 130: D 22 1 -> 2 (no-slip walls in the implicit momentum diffusion; D 22 1 is free slip for viscous flows)"<<endl;
        D22 = 2;
        }
    }
    
    // CFD: the legacy rheology source terms W 110 2 - 7 are retired (inconsistent units and directions, only with the
    // void inflow): the viscosity form W 110 1, or the rheology core W 130
    if(mpirank==0 && A10==6 && W90>0 && W110>=2)
    {
    cout<<"!!! ctrl.txt: W 110 "<<W110<<" is retired: the rheology acts through the viscosity (W 110 1), see the rheology core W 130 !!!"<<endl;
    exit(1);
    }
    
    // NHFLOW: A 519 not given in ctrl.txt -> bed friction (log law in the bottom layer) with a turbulence model (A 560)
    // or sediment transport (S 10): without it the flow slips over the bed, the velocity profile of an open channel
    // stays a plug and the bed shear stress of the sediment is too high; no bed friction otherwise (the previous default)
    if(mpirank==0 && A519<0)
    A519 = (A560>0 || S10>0) ? 1 : 0;

#ifndef REEF3D_USE_HYPRE
    // built without hypre (the default; hypre is opt-in via make HYPRE=1): the hypre solvers N 10 10-39 are replaced by REEFMG (N 10 1)
    if(mpirank==0 && N10>=10)
    {
        std::cout<<"N 10 "<<N10<<": this build has no hypre (hypre is opt-in: make HYPRE=1), using REEFMG (N 10 1)"<<std::endl;
        N10=1;
    }
#endif

    Iarray(ictrl,ctrlsize);
    Darray(dctrl,ctrlsize);

    if(mpirank==0)
    control::ctrlsend();

    pgc->globalctrl(this);

    if(mpirank>0)
    control::ctrlrecv();

    del_Iarray(ictrl,ctrlsize);
    del_Darray(dctrl,ctrlsize);

    read_grid();

    control::parse(this);

    lexer_ini();
}
