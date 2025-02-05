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

#include"sediment_part.h"
#include"partres.h"
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

void sediment_part::sediment_logic(lexer *p, fdm *a,ghostcell *pgc, turbulence *pturb)
{
    pst = new partres(p,pgc);
    
    s = new sediment_fdm(p);
    
    if(p->S90==0)
    pslide=new sandslide_v(p);   
    
    if(p->S90==1)
    pslide=new sandslide_f(p);
    
    if(p->S90==2)
    pslide=new sandslide_f2(p);
    
    if(p->S90==3)
    pslide=new sandslide_f3(p);
    
    if(p->S90==4)
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
    
    ptopo = new sediment_exner(p,pgc);
   
    

	p->gcin4a_count=p->gcin_count;
	p->gcout4a_count=p->gcout_count;
	
    
    prelax = new topo_relax(p);
	
	pbedshear  = new bedshear(p,pturb);
    
    volume_token=0;
    
    
    
}