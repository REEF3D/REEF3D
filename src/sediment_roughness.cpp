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

#include"sediment_roughness.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

sediment_roughness::sediment_roughness(lexer *p)
{

}

sediment_roughness::~sediment_roughness()
{
}

void sediment_roughness::start(lexer* p, ghostcell* pgc, sediment_fdm *s, slice &WL)
{
    double Delta,lambda,Ti;
    
    if(p->S36==0)
    SLICELOOP4
    s->ks_eff(i,j) = p->S21*p->S20;
    
    
    if(p->S36==1)
    SLICELOOP4
    {
    Ti  = (s->tau_i(i,j)-s->tau_crit(i,j))/s->tau_crit(i,j);
    
    Delta = WL(i,j) * 0.11 * pow(p->S20/WL(i,j),0.3) * (1.0 - pow(EE,-0.5*Ti)) * (25.0 - Ti);
     
    lambda = 7.3 * WL(i,j);
    
    if(Ti < 25.0 && Ti>0.0)
    s->ks_eff(i,j) = p->S21*p->S20 + 1.1 * Delta*(1.0-pow(EE,-25.0*Delta/lambda));
    
     
     
    if(Ti >= 25.0 || Ti<=0.0)
    s->ks_eff(i,j) = p->S21*p->S20;
    }
}
