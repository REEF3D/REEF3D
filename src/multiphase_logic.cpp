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

#include"multiphase_f.h"
#include"lexer.h"
#include"freesurface_header.h"
#include"multiphase_fluid_update_f.h"
#include"multiphase_fluid_update_rheology.h"
#include"heat_void.h"
#include"concentration_void.h"
#include"print_wsf.h"
#include"convection_header.h"

void multiphase_f::logic(lexer *p, fdm *a, ghostcell *pgc)
{
    pwsf1 = new print_wsf(p,a,pgc,1);
    pwsf2 = new print_wsf(p,a,pgc,2);

    pheat = new heat_void();
    pconc = new concentration_void();
    
    // Free Surface
    if(p->F300==0)
    {
        pfsf1 = new levelset_void(p,a,pgc,pheat,pconc);
        pfsf2 = new levelset_void(p,a,pgc,pheat,pconc);
    }
    else if(p->F300==1)
    {
        pfsf1 = new levelset_AB2(p,a,pgc,pheat,pconc);
        pfsf2 = new levelset_AB2(p,a,pgc,pheat,pconc);
    }
    else if(p->F300==2)
    {
        pfsf1 = new levelset_RK2(p,a,pgc,pheat,pconc);
        pfsf2 = new levelset_RK2(p,a,pgc,pheat,pconc);
    }
    else if(p->F300==3)
    {
        pfsf1 = new levelset_RK3(p,a,pgc,pheat,pconc);
        pfsf2 = new levelset_RK3(p,a,pgc,pheat,pconc);
    }

    //  Convection LSM
    if(p->F305==0)
        pmpconvec = new convection_void(p);
    else if(p->F305==1)
        pmpconvec = new fou(p);
    else if(p->F305==2)
        pmpconvec = new cds2(p);
    else if(p->F305==3)
        pmpconvec = new quick(p);
    else if(p->F305==4)
        pmpconvec = new weno_flux_nug(p);
    else if(p->F305==5)
        pmpconvec = new weno_hj_nug(p);
    else if(p->F305==6)
        pmpconvec = new cds4(p);
    else if((p->F305>=10 && p->F305<30) || (p->F305>=40 && p->F305<50))
        pmpconvec = new hires(p,p->F305);
    
    if(p->F310==0)
        preini = new reini_void(p);
    else if(p->F310==3)
        preini = new reini_RK3(p,1);
    else if(p->F310==11 || p->F310==13 || p->F310==14)
        preini = new directreini(p,a);
    
    if(p->F31==0)
        ppls = new particle_pls_void();
    else if(p->F31==1 || p->F31==2)
        ppls = new particle_pls(p,a,pgc);
    
    if(p->W90==0)
        pupdate = new multiphase_fluid_update_f(p,a,pgc);
    else if(p->W90>0)
        pupdate = new multiphase_fluid_update_rheology(p);
}
