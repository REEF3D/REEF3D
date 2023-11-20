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
    
    alpha[0] = 8.0/15.0;
    alpha[1] = 2.0/15.0;
    alpha[2] = 2.0/6.0;
    
    if(p->N40==23 || p->N40==33)
    {
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
    }
    
    if(p->N40==2 || p->N40==22)
    {
    alpha[0] = 1.0;
    alpha[1] = 0.5;
    }
    
    gamma[0] = 8.0/15.0;
    gamma[1] = 5.0/12.0;
    gamma[2] = 3.0/4.0;
    
    zeta[0] = 0.0;
    zeta[1] = -17.0/60.0;
    zeta[2] = -5.0/12.0;
}

sixdof_df_object::~sixdof_df_object(){}
    
