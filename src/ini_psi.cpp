/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::inipsi(lexer* p, fdm *a, ghostcell* pgc)
{
    double psim;
    int count;
    
    if(p->j_dir==0)        
    p->psi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    p->psi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
    
    
    if(p->B90>0 || p->B60>0)
    {
    // psi
        count=0;
        psim=0.0;
        LOOP
        if(fabs(a->phi(i,j,k))<5.0*p->DTM)
        {
        psim += p->DZN[KP];
        ++count;
        }
        
        count=pgc->globalisum(count);
        psim=pgc->globalsum(psim);
        
        p->psi = p->F45*psim/double(count);

    }
    
    p->psi0=p->psi;
    
    //cout<<p->mpirank<<" PSI0: "<<p->psi0<<" PSI: "<<p->psi<<" DTM: "<<p->DTM<<endl;
    
}