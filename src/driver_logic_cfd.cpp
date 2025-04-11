/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"reinisolid_RK3.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"FSI_header.h"
#include"vrans_header.h"

void driver::logic_cfd()
{
	makegrid_cds();
	pini = new initialize(p);

    if(p->mpirank==0)
        cout<<"starting ini"<<endl;
    pini->start(a,p,pgc);

    if(p->mpirank==0)
        cout<<"creating objects"<<endl;

    // time stepping
    if(p->N48==0)
        ptstep = new fixtimestep(p);
    else if(p->N48==1  && (p->D20!=0 && p->D20!=2))
        ptstep = new etimestep(p);
    else if(p->N48==1 && (p->D20==0 || p->D20>=2))
        ptstep = new ietimestep(p);
    
    assign_flux();
    // Multiphase
    if(p->F300==0)
        pmp = new multiphase_v();
    else if(p->F300>0)
        pmp = new multiphase_f(p,a,pgc,pflux);

    //discretization schemes
    //Convection
    if(p->D10==0)
        pconvec = new convection_void(p);
    else if(p->D10==1)
        pconvec = new fou(pflux);
    else if(p->D10==2)
        pconvec = new cds2(pflux);
    else if(p->D10==60)
        pconvec = new hcds6(pflux);
    else if(p->D10==3)
        pconvec = new quick(pflux);
    else if(p->D10==4)
        pconvec = new weno_flux_nug(p,pflux);
    else if(p->D10==5)
        pconvec = new weno_hj_nug(p);
    else if(p->D10==6)
        pconvec = new cds4(pflux);
    else if(p->D10==7)
        pconvec = new weno3_flux(p,pflux);
    else if(p->D10==8)
        pconvec = new weno3_hj(p);
    else if(p->D10>=10 && p->D10<30)
        pconvec = new hires(p,pflux,p->D10);

    // Convection Turbulence
    if(p->T12==0)
        pturbdisc = new convection_void(p);
    else if(p->T12==1)
        pturbdisc = new ifou(pflux);
    else if(p->T12==5)
        pturbdisc = new iweno_hj_df_nug(p);
    else if(p->T12==55)
        pturbdisc = new iweno_hj(p);

    //  Convection FSF && VOF
    if(p->F35==0 && p->F85==0)
        pfsfdisc = new convection_void(p);
    else if(p->F35==1 || p->F85==1)
        pfsfdisc = new fou(pflux);
    else if(p->F35==2 || p->F85==2)
        pfsfdisc = new cds2_alt(p);
    else if(p->F35==3 || p->F85==3)
        pfsfdisc = new quick(pflux);
    else if(p->F35==4 || p->F85==4)
        pfsfdisc = new weno_flux_nug(p,pflux);
    else if(p->F35==5)
        pfsfdisc = new weno_hj_df_nug(p);
    else if(p->F85==5)
        pfsfdisc = new weno_hj_nug(p);
    else if(p->F35==6 || p->F85==6)
        pfsfdisc = new cds4(pflux);
    else if(p->F35==7)
        pfsfdisc = new weno3_flux(p,pflux);
    else if(p->F35==8)
        pfsfdisc = new weno3_hj(p);
    else if((p->F35>=10 && p->F35<30) || (p->F35>=40 && p->F35<50))
        pfsfdisc = new hires(p,pflux,p->F35);
    else if((p->F85>=10 && p->F85<30) || (p->F85>=40 && p->F85<50))
        pfsfdisc = new hires(p,pflux,p->F85);
    else if(p->F85==51)
        pfsfdisc = new hric(pflux);
    else if(p->F85==52)
        pfsfdisc = new hric_mod(pflux);
    else if(p->F85==53)
        pfsfdisc = new cicsam(pflux);

    //  Convection Concentration
    if(p->C15==0)
        pconcdisc = new convection_void(p);
    else if(p->C15==1)
        pconcdisc = new fou(pflux);
    else if(p->C15==2)
        pconcdisc = new cds2_alt(p);
    else if(p->C15==3)
        pconcdisc = new quick(pflux);
    else if(p->C15==4)
        pconcdisc = new weno_flux_nug(p,pflux);
    else if(p->C15==5)
        pconcdisc = new weno_hj_nug(p);
    else if(p->C15==6)
        pconcdisc = new cds4(pflux);
    else if(p->C15==7)
        pconcdisc = new weno3_flux(p,pflux);
    else if(p->C15==8)
        pconcdisc = new weno3_hj(p);
    else if((p->C15>=10 && p->C15<30) || (p->C15>=40 && p->C15<50))
        pconcdisc = new hires(p,pflux,p->C15);
    else if(p->S60>0 && p->S60<10)
        pconcdisc = new weno_hj(p);
    else if(p->S60==11 || p->S60==12)
        pconcdisc = new iweno_hj(p);
    
    //turbulence model
    if(p->T10==0)
        pturb = new kepsilon_void();
    //ke
    else if(p->T10==1 || p->T10==21)
        pturb = new kepsilon_IM1(p,a,pgc);
    //kw
    else if(p->T10==2 || p->T10==22)
        pturb = new komega_IM1(p,a,pgc);
    //EARSM
    else if(p->T10==12)
        pturb = new EARSM_kw_IM1(p,a,pgc);
    // LES
    else if(p->T10==31)
        pturb = new LES_smagorinsky(p,a);
    else if(p->T10==33)
        pturb = new LES_WALE(p,a);

    //Heat
    if(p->H10==0)
        pheat = new heat_void();
    else if(p->H10==1)
        pheat = new heat_AB(p,a,pgc,pheat);
    else if(p->H10==2)
        pheat = new heat_RK2(p,a,pgc,pheat);
    else if(p->H10==3)
        pheat = new heat_RK3(p,a,pgc,pheat);
    else if(p->H10==4)
        pheat = new heat_RK3CN(p,a,pgc,pheat);

    //Convection Heat
    if(p->H15==0)
        pheatdisc = new convection_void(p);
    else if(p->H15==1)
        pheatdisc = new fou(pflux);
    else if(p->H15==2)
        pheatdisc = new cds2(pflux);
    else if(p->H15==60)
        pheatdisc = new hcds6(pflux);
    else if(p->H15==3)
        pheatdisc = new quick(pflux);
    else if(p->H15==4)
        pheatdisc = new weno_flux_nug(p,pflux);
    else if(p->H15==5)
        pheatdisc = new weno_hj_nug(p);
    else if(p->H15==6)
        pheatdisc = new cds4(pflux);
    else if(p->H15==7)
        pheatdisc = new weno3_flux(p,pflux);
    else if(p->H15==8)
        pheatdisc = new weno3_hj(p);
    else if(p->H15==9)
        pheatdisc = new weno_flux(pflux);
    else if(p->H15>=10 && p->H15<30)
        pheatdisc = new hires(p,pflux,p->H15);

    // Concentration
    if(p->C10==0)
        pconc = new concentration_void();
    else if(p->C10==1)
        pconc = new concentration_AB(p,a,pgc);
    else if(p->C10==2)
        pconc = new concentration_RK2(p,a,pgc);
    else if(p->C10==3)
        pconc = new concentration_RK3(p,a,pgc);

    //Diffusion
    // momentum and scalars
    if(p->D20==0)
        pdiff = new diff_void();
    else if(p->D20==1)
        pdiff = new ediff2(p);
    else if(p->D20==2 && p->j_dir==0)
        pdiff = new idiff2_FS_2D(p);
    else if(p->D20==2 && p->j_dir==1)
        pdiff = new idiff2_FS(p);
    else if(p->D20==3 && p->j_dir==1)
        pdiff = new idiff2_CN(p);

    // turbulence
    if(p->D20==0 || p->T10==0)
        pturbdiff = new diff_void();
    else if(p->T10>0 && p->D20>0)
        pturbdiff = new idiff2(p);

    // concentration
    if(p->D20==0 || p->C10==0)
        pconcdiff = new diff_void();
    else if(p->D20==1 && p->C10<=10 && p->C10>0)
        pconcdiff = new ediff2(p);
    else if(p->D20>=2 && p->C10<=10 && p->C10>0)
        pconcdiff = new idiff2_FS(p);

    // Free Surface
    if((p->F30==0 && p->F80==0) || (p->N40==22||p->N40==23||p->N40==33))
        pfsf = new levelset_void(p,a,pgc,pheat,pconc);
    else if(p->F30==1 && p->N40!=22 && p->N40!=23 && p->N40!=33)
        pfsf = new levelset_AB2(p,a,pgc,pheat,pconc);
    else if(p->F30==2 && p->N40!=22 && p->N40!=23 && p->N40!=33)
        pfsf = new levelset_RK2(p,a,pgc,pheat,pconc);
    else if(p->F30==3 && p->N40!=22 && p->N40!=23 && p->N40!=33)
        pfsf = new levelset_RK3(p,a,pgc,pheat,pconc);
    else if(p->F80==1)
        pfsf = new VOF_AB(p,a,pgc,pheat,pflux);
    else if(p->F80==3)
        pfsf = new VOF_RK3(p,a,pgc,pheat,pflux);
    else if(p->F80==4)
        pfsf = new VOF_PLIC(p,a,pgc,pheat);

    // Reini
    if(p->F40==0)
        preini = new reini_void(p);
    else if(p->F40==3 || p->F40==23)
        preini = new reini_RK3(p,1);
    else if(p->F40==11)
        preini = new directreini(p,a);

    // Particle based level set
    if(p->F31==0)
        ppls = new particle_pls_void();
    else if(p->F31==1 || p->F31==2)
        ppls = new particle_pls(p,a,pgc);

    assign_density();

    // Pressure scheme
    if(p->D30==0)
        ppress = new pressure_void();
    else if((p->D30==1 || p->D30==2 || p->D30==3))
        ppress = new pjm_corr(p,pd);
    else if(p->D30==10)
        ppress = new pjm_hydrostatic(pd);

    //Solver
    assign_solver();

    //Poison Solver
    assign_poisson_solver();

    //VRANS
    assign_VRANS();

    //IOFlow
    assign_IOFlow();

    //Potential Flow Solver
    if(p->I11==0)
        potflow = new potential_v();
    else if(p->I11==1)
        potflow = new potential_f(p);
    else if(p->I11==2)
        potflow = new potential_water(p);

    // Benchmark
    if(p->F150==0)
        pbench = new benchmark_void();
    else if(p->F150==1)
        pbench = new benchmark_vortex(p,a);
    else if(p->F150==2)
        pbench = new benchmark_disk(p,a);
    else if(p->F150==3)
        pbench = new benchmark_vortex3D(p,a);
    else if(p->F150==4)
        pbench = new benchmark_TaylorGreen(p,a);
    else if(p->F150==11)
        pbench = new benchmark_convection(p,a);

    // Printer
    if(p->P10==2)
        pprint = new vtr3D(p,a,pgc);
    else if(p->P10==3)
        pprint = new vts3D(p,a,pgc);
    else
        pprint = new vtu3D(p,a,pgc);

    // Data
    assign_data();

    // Sediment
    if(p->S10>0)
    {
        if(p->Q10==0)
            psed = new sediment_f(p,pgc,pturb,pBC);
        if(p->Q10==1)
            psed = new sediment_part(p,a,pgc,pturb,pBC);
    }
    else
        psed = new sediment_void();

    // Reinitopo
    if(p->S10>0 || p->G1==1 || p->toporead==1)
        assign_reinitopo();

    // Reinisolid
    if(p->solidread==0 || p->G40==0)
        preso = new reinitopo_void();
    else if(p->solidread==1 && p->G40>0)
        preso = new reinisolid_RK3(p);
    
    // 6DOF
    if(p->X10==0)
        p6dof = new sixdof_void(p,pgc);
    else if(p->X10==1)
        p6dof = new sixdof_cfd(p,a,pgc);

    // FSI
    if(p->Z10==0)
        pfsi = new fsi_void(p,pgc);
    else if(p->Z10==1)
        pfsi = new fsi_strips(p,pgc);

    // Velocities
    if(p->N40==0)
        pmom = new momentum_void();
    else if(p->N40==2)
        pmom = new momentum_RK2(p,a,pconvec,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow,pfsi);
    else if(p->N40==3)
        pmom = new momentum_RK3(p,a,pconvec,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow,pfsi);
    else if(p->N40==4)
    {
        pmom = new momentum_void();
        if(p->X10==0 && p->Z10==0)
            pmom_sf = new momentum_RKLS3_sf(p,a,pgc,pconvec,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow);
        else if(p->X10==1 || p->Z10>0)
            pmom_df = new momentum_RKLS3_df(p,a,pgc,pconvec,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow);
    }
    else if(p->N40==5)
        pmom = new momentum_RK3CN(p,a,pconvec,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow,pfsi);
    else if(p->N40==22)
        pmom = new momentum_FC2(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
    else if(p->N40==23)
        pmom = new momentum_FC3(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
    else if(p->N40==33)
        pmom = new momentum_FCC3(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi,pd);
}

void driver::patchBC_logic()
{
    if((p->B440>0 || p->B441>0 || p->B442>0) && p->A10==2)
        pBC = new patchBC_2D(p,pgc);
    else if((p->B440>0 || p->B441>0 || p->B442>0) && p->A10==6)
        pBC = new patchBC(p,pgc);
    else if((p->B440==0 && p->B441==0 && p->B442==0) || (p->A10!=2 && p->A10!=6))
        pBC = new patchBC_void(p);
}
