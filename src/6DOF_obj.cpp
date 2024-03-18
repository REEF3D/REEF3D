/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"
#include"6DOF_motionext_fixed.h"
#include"6DOF_motionext_file.h"
#include"6DOF_motionext_void.h"

sixdof_obj::sixdof_obj(lexer *p, ghostcell *pgc, int number) : gradient(p), dt(p), L(p), 
                                                                                f(p), frk1(p), cutl(p), cutr(p), 
                                                                                fbio(p),n6DOF(number),
                                                                                epsifb(1.6*p->DXM), epsi(1.6),vertice(p),
                                                                                nodeflag(p),interfac(1.6),zero(0.0),eta(p)
{
    prdisc = new reinidisc_fsf(p);
    
    triangle_token=0;
    printnormal_count=0;
    
    alpha[0] = 8.0/15.0;
    alpha[1] = 2.0/15.0;
    alpha[2] = 2.0/6.0;
    
    gamma[0] = 8.0/15.0;
    gamma[1] = 5.0/12.0;
    gamma[2] = 3.0/4.0;
    
    zeta[0] = 0.0;
    zeta[1] = -17.0/60.0;
    zeta[2] = -5.0/12.0;
    
    if(p->N40==3 || p->N40==23 || p->N40==33)
    {
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
    
    gamma[0] = 0.0;
    gamma[1] = 0.0;
    gamma[2] = 0.0;
    
    zeta[0] = 0.0;
    zeta[1] = 0.0;
    zeta[2] = 0.0;
    }
    
    if(p->N40==2 || p->N40==22)
    {
    alpha[0] = 1.0;
    alpha[1] = 0.5;
    
    gamma[0] = 0.0;
    gamma[1] = 0.0;
    gamma[2] = 0.0;
    
    zeta[0] = 0.0;
    zeta[1] = 0.0;
    zeta[2] = 0.0;
    }
    
    
    if(p->X210==0 && p->X211==0)
    pmotion = new sixdof_motionext_void(p,pgc);
    
    if(p->X210==1 || p->X211==1)
    pmotion = new sixdof_motionext_fixed(p,pgc);
    
    if(p->X240>0)
    pmotion = new sixdof_motionext_file(p,pgc);
}

sixdof_obj::~sixdof_obj()
{
}
    
