/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"

sixdof_df_object::sixdof_df_object(lexer *p, fdm *a, ghostcell *pgc,int number) : gradient(p), dt(p), L(p), 
                                                                                f(p), frk1(p), cutl(p), cutr(p), 
                                                                                fbio(p),n6DOF(number),
                                                                                epsifb(1.6*p->DXM), epsi(1.6)
{
    prdisc = new reinidisc_fsf(p);
    
    triangle_token=0;
    printnormal_count=0;
}

sixdof_df_object::~sixdof_df_object(){}
    
