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

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"sediment_exner.h"
#include"reinitopo.h"
#include"suspended.h"
#include"bedload_VR.h"
#include"bedload_einstein.h"
#include"bedload_MPM.h"
#include"bedload_EF.h"
#include"bedload_EH.h"
#include"bedload_void.h"
#include"bedconc_void.h"
#include"bedconc_VR.h"
#include"bedshear.h"
#include"sandslide_f.h"
#include"sandslide_f2.h"
#include"sandslide_pde.h"
#include"sandslide_steepest_descent.h"
#include"sandslide_weighted_multidir.h"
#include"sandslide_v.h"
#include"topo_relax.h"
#include"vrans_v.h"
#include"vrans_f.h"
#include"bedslope.h"
#include"reduction_void.h"
#include"reduction_parker.h"
#include"reduction_deyemp.h"
#include"reduction_deyana.h"
#include"reduction_FD.h"
#include"reduction_FD_gamma.h"
#include"diff_void.h"
#include"ediff2.h"
#include"idiff2.h"
#include"idiff2_FS.h"
#include"idiff2_FS_2D.h"
#include"convection_void.h"
#include"weno_hj_nug.h"
#include"ifou.h"
#include"suspended_void.h"
#include"suspended_RK2.h"
#include"suspended_RK3.h"
#include"suspended_IM1.h"
#include"bedload_direction_f.h"
#include"bedload_direction_v.h"
#include"nhflow_diff_void.h"
#include"nhflow_ediff.h"
#include"nhflow_idiff.h"
#include"nhflow_idiff_2D.h"
#include"nhflow_scalar_iweno.h"
#include"nhflow_scalar_ifou.h"
#include"nhflow_scalar_void.h"
#include"nhflow_suspended_void.h"
#include"nhflow_suspended_IM1.h"
#include"bedprobe_point.h"
#include"bedprobe_max.h"
#include"bedshear_probe.h"
#include"bedshear_max.h"
#include"bedprobe_line_x.h"
#include"bedprobe_line_y.h"
#include"sediment_roughness.h"
#include"sediment_mixture.h"
#include"sediment_waves.h"
#include"bedload_waves.h"
#include"sediment_susp2D.h"
#include"slice1.h"
#include"slice2.h"

void sediment_f::sediment_logic(lexer *p, ghostcell *pgc, turbulence *pturb)
{
    // bend morphology (sediment_exner_bend.cpp): the spiral-flow correction is for the depth-averaged
    // SFLOW (NHFLOW and CFD resolve the secondary flow); both corrections act on the conservative Exner forms
    if((p->S95>0 && p->A10!=2) || ((p->S86>0 || p->S95>0) && p->S31==3))
    {
        if(p->mpirank==0)
        cout<<"!!! S 95 (spiral-flow correction) is for SFLOW (A 10 2) only; S 86 and S 95 need S 31 1 or 2 !!!"<<endl;
        exit(1);
    }
    
    // the lags of the bend morphology (sediment_exner_lag.cpp): S 95 0, 1 or 2; the suspended part S 64 of the
    // depth-averaged transport for SFLOW with one fraction, 0 <= S 64 <= 1
    if(p->S95<0 || p->S95>2 || p->S97<0.0
       || (p->S64!=0.0 && (p->A10!=2 || p->S51>0 || p->S64<0.0 || p->S64>1.0 || p->S65<=0.0 || p->S31==3)))
    {
        if(p->mpirank==0)
        cout<<"!!! S 95 is 0, 1 or 2 (S 97 >= 0); S 64 (suspended part, 0 .. 1, S 65 > 0) is for SFLOW (A 10 2) with one fraction and S 31 1 or 2 !!!"<<endl;
        exit(1);
    }

    s = new sediment_fdm(p);
    
    // multi-fraction bed
    if(p->S51>0)
    s->pmix = new sediment_mixture(p);
    
    
    if(p->S11==0)
    pbed = new bedload_void();

    if(p->S11==1)
    pbed = new bedload_VR(p);

    if(p->S11==2)
    pbed = new bedload_MPM(p);
	
	if(p->S11==3)
    pbed = new bedload_EF(p);
    
    if(p->S11==4)
    pbed = new bedload_EH(p);
    
    if(p->S11==5)
    pbed = new bedload_einstein(p);
    
    // SFLOW wave formulas (bedload_waves.h)
    if(p->S11>=6 && p->S11<=9)
    {
    check_sflow_waves(p);
    pbedw = new bedload_waves(p);
    pbed = pbedw;
    
    // the flow velocities of the sediment step and the transport velocity (transport_direction)
    Psave = new slice1(p);
    Qsave = new slice2(p);
    Vtx = new slice4(p);
    Vty = new slice4(p);
    }
    
    if(p->S12==0)
    pcbed = new bedconc_void(p);
    
    if(p->S12==1)
    pcbed = new bedconc_VR(p);
    
    
    if(p->S90==0)
    pslide=new sandslide_v(p);   
    
    if(p->S90==1)
    pslide=new sandslide_f(p);
    
    if(p->S90==2)
    pslide=new sandslide_f2(p);
    
    if(p->S90==3)
    pslide=new sandslide_steepest_descent(p);
    
    if(p->S90==4)
    pslide=new sandslide_weighted_multidir(p);
    
    if(p->S90==5)
    pslide=new sandslide_pde(p);
    
    
    
    if(p->S10!=2 && p->A10==6)
	pvrans = new vrans_v(p,pgc);
	
	if(p->S10==2 && p->A10==6)
	pvrans = new vrans_f(p,pgc);
    
    
    pslope = new bedslope(p);
    
    if(p->S80==0)
    preduce=new reduction_void(p);

    if(p->S80==1)
    preduce=new reduction_parker(p);

    if(p->S80==2)
    preduce=new reduction_deyemp(p);

    if(p->S80==3)
    preduce=new reduction_deyana(p);
	
	if(p->S80==4)
    preduce=new reduction_FD(p);
    
    if(p->S80==5)
    preduce=new reduction_FD_gamma(p);
    
    ptopo = new sediment_exner(p,pgc);
    
    // Exner formulations that do not conserve the bed volume
    if(p->mpirank==0 && amrgrid==0)
    {
        if(p->S31==3)
        cout<<"WARNING S 31 3: advective Exner form sign(u)*dqb/dx, not conservative where the transport direction turns; use S 31 1 or 2"<<endl;
        
        if(p->S32==3 && p->S31==1)
        cout<<"WARNING S 32 3 with S 31 1: advective (HJ) form, not conservative; use S 32 1, 2 or 4"<<endl;
        
        if(p->S32==5)
        cout<<"WARNING S 32 5: WENO-HJ derivative, not conservative; use S 32 1, 2 or 4"<<endl;
        
        if(p->S61==2)
        cout<<"NOTE S 61 2: the reference concentration is evaluated at the first cell centre (bedconc_VR), the near-bed concentration is taken there as well (S 61 1)"<<endl;
        
        if(p->S62==2 && p->A10!=5)
        cout<<"WARNING S 62 2: the suspended flux qbs is only computed by NHFLOW, S 62 1 (exchange E - D in the Exner equation) is used"<<endl;
        
        if(p->S62==2 && p->A10==5 && p->S12>0)
        cout<<"NOTE S 62 2: the bed sees div(qbs) while the suspended solver still exchanges E - D with the water column, suspended sediment is counted twice; S 62 1 is the conservative coupling"<<endl;
        
        if(p->S80>0 && p->S85==1)
        cout<<"NOTE S 80 "<<p->S80<<" with S 85 1: the bed slope reduces the critical shear stress and scales the bedload (1-1.3 dz/ds), slope effect counted twice; S 85 0 uses the critical shear stress reduction only"<<endl;
    }

    // S 62 2 needs qbs, which only NHFLOW computes
    if(p->S62==2 && p->A10!=5)
    p->S62=1;
    
    // Suspended Sediments
    // Suspended NHFLOW
    if(p->A10==5)
    {  
    if(p->S12==0)
	pnhfsuspdiff=new nhflow_diff_void(p);
    
    if(p->S12==0)
    pnhfsusp = new nhflow_suspended_void(p);
    
    if(p->S12==0)
	pnhfsuspdisc=new nhflow_scalar_void(p);
	
    
    if(p->S12>0 && p->j_dir==1)
    pnhfsuspdiff = new nhflow_idiff(p);
        
    if(p->S12>0 && p->j_dir==0)
    pnhfsuspdiff = new nhflow_idiff_2D(p);
    
	if(p->S12>0)
	pnhfsuspdisc= new nhflow_scalar_ifou(p,2);   // conservative form for D*C, WVEL = omegaF - ws (settling as a face flux)
    
    if(p->S12>0)
    pnhfsusp = new nhflow_suspended_IM1(p);
    }
    
    // Suspended CFD
    if(p->A10==6)
    {  
    if(p->S12==0)
	psuspdiff=new diff_void();
    
    if(p->S12==0)
    psusp = new suspended_void();
    
    if(p->S12==0)
	psuspdisc=new convection_void(p);
    
    
    if(p->S12>=1)
    psuspdiff=new idiff2(p);
    
    if(p->S12>=1)
	psuspdisc=new ifou(p);
    

    if(p->S12>=1)
    psusp = new suspended_IM1(p);
    }
    
    if(p->S85==0)
    pbeddir = new bedload_direction_v(p);
    
    if(p->S85==1)
    pbeddir = new bedload_direction_f(p);
    
	
	
    
    prelax = new topo_relax(p);
    
    pbedshear  = new bedshear(p,pturb);
	
    
    volume_token=0;
    
    
    // the grids of the mesh refinement and the blocks: no probes (level 0 writes them)
    if(amrgrid==0)
    {
    if(p->P121>0)
	pbedpt = new bedprobe_point(p,pgc,s);

	if(p->P122>0)
	pbedmax = new bedprobe_max(p,pgc,s);

	if(p->P123>0)
	pbedlinex=new bedprobe_line_x(p,pgc,s);

	if(p->P124>0)
	pbedliney=new bedprobe_line_y(p,pgc,s);

	if(p->P125>0)
	pbedshearprobe = new bedshear_probe(p,pgc);

	if(p->P126>0)
	pbedshearmax = new bedshear_max(p,pgc);
    }
    
    pks = new sediment_roughness(p);
    
    // SFLOW, NHFLOW: waves in the sediment transport (sediment_waves.h); SFLOW: depth-averaged suspended load (sediment_susp2D.h)
    if(p->S110>0 || p->S120>0)
    check_sflow_waves(p);
    
    if(p->S110>0)
    {
    pwaves = new sediment_waves(p);
    s->pwave = pwaves;
    }
    
    if(p->S120>0)
    psusp2D = new sediment_susp2D(p);
    
    // suspended load (SFLOW S 120, NHFLOW S 12): the bed takes S 35 dtsed/T of the exchange of the flow time T since the
    // last bed update
    if((p->S120>0 || (p->A10==5 && p->S12>0 && p->S62==1)) && p->mpirank==0 && amrgrid==0 && !(p->S15==4 || (p->S15==3 && p->S42==1 && p->S44==1)))
    cout<<"NOTE suspended load: the bed takes S 35 dtsed/T of the suspended exchange of the flow time T since the last bed update; "
        <<"bed and suspended load conserve the sediment with S 15 4 or S 15 3, S 42 1, S 44 1 (and S 17 1, S 35 1)"<<endl;
}

void sediment_f::check_sflow_waves(lexer *p)
{
    int stop=0;
    
    if(p->S120>0 && p->A10!=2)
    {
    if(p->mpirank==0)
    cout<<"!!! S 120 (depth-averaged suspended load) is for SFLOW (A 10 2), NHFLOW has the suspended load S 12 !!!"<<endl;
    stop=1;
    }
    
    if((p->S110>0 || (p->S11>=6 && p->S11<=9)) && p->A10!=2 && p->A10!=5)
    {
    if(p->mpirank==0)
    cout<<"!!! S 11 6-9 and S 110 are for SFLOW (A 10 2) and NHFLOW (A 10 5) !!!"<<endl;
    stop=1;
    }
    
    // NHFLOW: the formulas and wave parameters are averaged over the wave cycle, for NHFLOW as the current model of
    // phase-averaged waves (SEASTATE or prescribed); waves resolved by NHFLOW are already in its velocity
    if(p->A10==5 && (p->S110>0 || (p->S11>=6 && p->S11<=9)) && (p->B90>0 || p->B98>=2))
    {
    if(p->mpirank==0)
    cout<<"!!! S 11 6-9 and S 110 with NHFLOW: for phase-averaged waves (A 750 1 or S 110 2), not with the waves "
        <<"generated by NHFLOW (B 90, B 98), whose orbital motion is in the near-bed velocity already (use S 11 1-5, S 16 9) !!!"<<endl;
    stop=1;
    }
    
    if(p->A10==5 && p->S110>0 && p->S16==9)
    {
    if(p->mpirank==0)
    cout<<"!!! S 110 and S 16 9 both add the waves to the bed shear stress: use one of them !!!"<<endl;
    stop=1;
    }
    
    if(p->S11>=6 && p->S51>0)
    {
    if(p->mpirank==0)
    cout<<"!!! S 11 6-9: one sediment fraction (S 51 0) !!!"<<endl;
    stop=1;
    }
    
    if(p->S110==1 && p->A750!=1)
    {
    if(p->mpirank==0)
    cout<<"!!! S 110 1 needs REEF3D::SEASTATE coupled with SFLOW or NHFLOW (A 750 1) !!!"<<endl;
    stop=1;
    }
    
    if(p->S110==2 && (!(p->S111>0.0) || !(p->S112>0.0)))
    {
    if(p->mpirank==0)
    cout<<"!!! S 110 2 needs Hs (S 111) and Tp (S 112) > 0 !!!"<<endl;
    stop=1;
    }
    
    if(p->S110<0 || p->S110>2 || p->S120<0 || p->S120>2 || (p->S11>=10))
    {
    if(p->mpirank==0)
    cout<<"!!! S 11 0-9, S 110 0-2, S 120 0-2 !!!"<<endl;
    stop=1;
    }
    
    if(p->S120==2 && p->S11!=6 && p->S11!=8)
    {
    if(p->mpirank==0)
    cout<<"!!! S 120 2 takes C_eq from the suspended part of S 11 6 or 8 !!!"<<endl;
    stop=1;
    }
    
    if(p->S120>0 && p->S64>0.0)
    {
    if(p->mpirank==0)
    cout<<"!!! S 120 (suspended load model) and S 64 (suspended part of the bed load) exclude each other !!!"<<endl;
    stop=1;
    }
    
    if(p->S11>=6 && p->S95>0)
    {
    if(p->mpirank==0)
    cout<<"!!! S 11 6-9 set the transport direction, S 95 (spiral flow) is not combined with them !!!"<<endl;
    stop=1;
    }
    
    // v1: Cartesian grid, one grid
    const int amr_on = ((p->G1>0 && (p->G10>0 || p->G20>0.0 || p->G22>0 || p->G12>0)) || p->geo_mb>0) ? 1 : 0;
    
    if(amr_on==1 || p->G50>0 || p->geo_curv==1)
    {
    if(p->mpirank==0)
    cout<<"!!! S 11 6-9, S 110 and S 120: Cartesian grid without mesh refinement or blocks (G 1, G 50, G 54) !!!"<<endl;
    stop=1;
    }
    
    if(stop==1)
    exit(1);
}