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