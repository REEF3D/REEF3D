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
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"driver.h"
#include"lexer.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_df.h"
#include"density_heat.h"
#include"density_vof.h"
#include"solver_header.h"
#include"ioflow_void.h"
#include"ioflow_f.h"
#include"iowave.h"
#include"ioflow_gravity.h"
#include"vrans_header.h"
#include"data_f.h"
#include"data_void.h"
#include"reinitopo_AB2.h"
#include"reinitopo_RK3.h"
#include"reinitopo_void.h"

void driver::assign_density()
{
    if(p->C10>0)
        pd = new density_conc(p,pconc);
    else if(p->F80>0)
        pd = new density_vof(p);
    else if(p->H10>0)
        pd = new density_heat(p,pheat);
    else if (p->W30==1)
        pd = new density_comp(p);
    else
        pd = new density_df(p);
}

void driver::assign_poisson_solver()
{
    if(p->N10==0)
        ppoissonsolv = new solver_void();
    else if(p->N10==1 && p->j_dir==0)
        ppoissonsolv = new bicgstab_ijk_2D(p);
    else if(p->N10==1 && p->j_dir==1)
        ppoissonsolv = new bicgstab_ijk(p);
    #ifdef HYPRE_COMPILATION
    else if(p->N10>=10 && p->N10<20)
        ppoissonsolv = new hypre_struct(p,pgc,p->N10,p->N11);
    else if(p->N10>=20 && p->N10<30)
        ppoissonsolv = new hypre_aij(p);
    else if(p->N10>=30 && p->N10<40)
        ppoissonsolv = new hypre_sstruct(p,pgc);
    #endif
}

void driver::assign_IOFlow()
{
    if(p->B60==0 && p->B90==0 && p->B180==0)
        pflow = new ioflow_v(p,pgc,pBC);
    else if(p->B60>=1)
        pflow = new ioflow_f(p,pgc,pBC);
    else if(p->B90>=1)
        pflow = new iowave(p,pgc,pBC);
    else if(p->B180==1||p->B191==1||p->B192==1)
        pflow = new ioflow_gravity(p,pgc,pBC);
}

void driver::assign_VRANS()
{
    if(p->B269==0)
        pvrans = new vrans_v(p,pgc);
    else if(p->B269==1)
        pvrans = new vrans_f(p,pgc);
    else if(p->B269==2)
        pvrans = new vrans_veg(p,pgc);
    else if(p->B269==3)
        pvrans = new vrans_net(p,pgc);
}

void driver::assign_solver()
{
    if(p->j_dir==0)
        psolv = new bicgstab_ijk_2D(p);
    else if(p->j_dir==1)
        psolv = new bicgstab_ijk(p);
}

void driver::assign_data()
{
    if(p->P150==0)
        pdata = new data_void();
    else if(p->P150>0)
        pdata = new data_f(p);
}

void driver::assign_reinitopo()
{
    if(p->G40==0)
        preto = new reinitopo_void();
    else if(p->G40==1)
        preto = new reinitopo_AB2(p);
    else if(p->G40==3)
        preto = new reinitopo_RK3(p);
}