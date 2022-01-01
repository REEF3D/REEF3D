/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
--------------------------------------------------------------------*/

#include"density_fsm.h"
#include"lexer.h"
#include"fdm.h"

density_fsm::density_fsm(lexer* p) : epsi(p->F45*p->DXM), eps(2.1*p->DXM)
{
        if(p->j_dir==0)        
        psi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
        if(p->j_dir==1)
        psi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
}

density_fsm::~density_fsm()
{
}

double density_fsm::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

    if(phival>psi)
    H=1.0;

    if(phival<-psi)
    H=0.0;

    if(fabs(phival)<=psi)
    H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
    
    if (aa == 1)
    {
        H_fb = a->fbh1(i,j,k);
    }
    else if (aa == -1)
    {
        H_fb = a->fbh1(i-1,j,k);
    }
    else if (bb == 1)
    {
        H_fb = a->fbh2(i,j,k);
    }
    else if (bb == -1)
    {
        H_fb = a->fbh2(i,j-1,k);
    }
    else if (cc == 1)
    {
        H_fb = a->fbh3(i,j,k);
    }
    else if (cc == -1)
    {
        H_fb = a->fbh3(i,j,k-1);
    }
    else
    {
        H_fb = a->fbh4(i,j,k);
    }
     
    roval = p->W_fb*H_fb + (1.0 - H_fb)*(p->W1*H + p->W3*(1.0 - H));

	return roval;		
}




