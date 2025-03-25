/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sediment_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"sediment_exner.h"
#include"reinitopo.h"
#include"suspended.h"
#include"bedload_VR.h"
#include"bedload_einstein.h"
#include"bedload_MPM.h"
#include"bedload_MPM.h"
#include"bedload_EF.h"
#include"bedload_EH.h"
#include"bedload_void.h"
#include"bedconc_void.h"
#include"bedconc_VR.h"
#include"bedshear.h"
#include"sandslide_f.h"
#include"sandslide_f2.h"
#include"sandslide_f3.h"
#include"sandslide_nz.h"
#include"sandslide_pde.h"
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
#include"iweno_hj_nug.h"
#include"suspended_void.h"
#include"suspended_RK2.h"
#include"suspended_RK3.h"
#include"suspended_IM1.h"
#include"bedload_direction_f.h"
#include"bedload_direction_v.h"

void sediment_f::sediment_logic(lexer *p, ghostcell *pgc, turbulence *pturb)
{
    s = new sediment_fdm(p);
    
    // Bedload
    switch(p->S11)
    {
        case 0:
        default:
            pbed = new bedload_void();
            break;
        case 1:
            pbed = new bedload_VR(p);
            break;
        case 2:
            pbed = new bedload_MPM(p);
            break;
        case 3:
            pbed = new bedload_EF(p);
            break;
        case 4:
            pbed = new bedload_EH(p);
            break;
        case 5:
            pbed = new bedload_einstein(p);
            break;
    }
    
    // Bed Concentration
    switch(p->S12)
    {
        case 0:
        default:
            pcbed = new bedconc_void(p);
            break;
        case 1:
            pcbed = new bedconc_VR(p);
            break;
    }
    
    // Sandslide
    switch(p->S90)
    {
        case 0:
        default:
            pslide = new sandslide_v(p);
            break;
        case 1:
            pslide = new sandslide_f(p);
            break;
        case 2:
            pslide = new sandslide_f2(p);
            break;
        case 3:
            pslide = new sandslide_f3(p);
            break;
        case 4:
            pslide = new sandslide_pde(p);
            break;
        case 5:
            pslide = new sandslide_nz(p);
            break;
    }
    
    // VRANS
    if(p->A10==6)
    {
        if(p->S10==2)
            pvrans = new vrans_f(p,pgc);
        else
            pvrans = new vrans_v(p,pgc);
    }
    else
        pvrans = nullptr;
    
    // Bed Slope
    pslope = new bedslope(p);
    
    // Reduction
    switch(p->S80)
    {
        case 0:
        default:
            preduce = new reduction_void(p);
            break;
        case 1:
            preduce = new reduction_parker(p);
            break;
        case 2:
            preduce = new reduction_deyemp(p);
            break;
        case 3:
            preduce = new reduction_deyana(p);
            break;
        case 4:
            preduce = new reduction_FD(p);
            break;
        case 5:
            preduce = new reduction_FD_gamma(p);
            break;
    }
    
    // Exner
    ptopo = new sediment_exner(p,pgc);
    
    // Suspended load
    if(p->A10!=6 || p->S60==0)
    {
        psuspdiff = new diff_void();
        psuspdisc = new convection_void(p);
        psusp = new suspended_void();
    }
    else if(p->A10==6 && p->S60>0)
    {
        // suspended diffusion
        if(p->S60<11 && p->j_dir==0)
            psuspdiff = new idiff2_FS_2D(p);
        else if(p->S60<11  && p->j_dir==1)
            psuspdiff = new idiff2_FS(p);
        else if(p->S60>10)
            psuspdiff = new idiff2(p);
        else
            pbeddir = nullptr;
        
        // suspended convection
        if(p->S60<11)
            psuspdisc = new weno_hj_nug(p);
        else if(p->S60>10)
            psuspdisc = new iweno_hj_nug(p);
        else
            pbeddir = nullptr;
        
        // suspended load
        if(p->S60==2)
            psusp = new suspended_RK2(p);
        else if(p->S60==3)
            psusp = new suspended_RK3(p);
        else if(p->S60==11)
            psusp = new suspended_IM1(p);
        else
            pbeddir = nullptr;
    }
    
    // Bedload direction    
    if(p->S85==1)
        pbeddir = new bedload_direction_f(p);
    else if(p->S85==0)
        pbeddir = new bedload_direction_v(p);
    else
        pbeddir = nullptr;
    
    // Topo Relax
    prelax = new topo_relax(p);
    
    // Bed Shear
    pbedshear = new bedshear(p,pturb);
    
    p->gcin4a_count=p->gcin_count;
    p->gcout4a_count=p->gcout_count;

    volume_token=0;
}
