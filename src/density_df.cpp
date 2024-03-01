/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"density_df.h"
#include"lexer.h"
#include"fdm.h"

density_df::density_df(lexer* p) : epsi(p->F45*p->DXM), eps(2.1*p->DXM)
{
    H=0.0;
}

density_df::~density_df()
{
}

double density_df::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double factor = 1.0;
    
    if(p->j_dir==0 && p->X46==1) 
    if(0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc)) <- 0.5*(1.0/2.0)*(p->DRM+p->DTM))
    factor = 2.0;
    
    if(p->j_dir==1 && p->X46==1) 
    if(0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc)) <- 0.5*(1.0/3.0)*(p->DRM+p->DSM+p->DTM))
    factor = 2.0;
    
    phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

    if(phival>factor*p->psi)
    H=1.0;

    if(phival<-factor*p->psi)
    H=0.0;

    if(fabs(phival)<=factor*p->psi)
    H=0.5*(1.0 + phival/(factor*p->psi) + (1.0/PI)*sin((PI*phival)/(factor*p->psi)));
    
    if(p->X15==1)
    {    
        if (aa == 1)
        H_fb = a->fbh1(i,j,k);
        
        else if (aa == -1)
        H_fb = a->fbh1(i-1,j,k);

        else if (bb == 1)
        H_fb = a->fbh2(i,j,k);

        else if (bb == -1)
        H_fb = a->fbh2(i,j-1,k);

        else if (cc == 1)
        H_fb = a->fbh3(i,j,k);

        else if (cc == -1)
        H_fb = a->fbh3(i,j,k-1);

        else
        H_fb = a->fbh4(i,j,k);
     
    roval = p->W_fb*H_fb + (1.0 - H_fb)*(p->W1*H + p->W3*(1.0 - H));
    }
    
    if(p->X15==2)
    roval = p->W1*H + p->W3*(1.0-H);

	return roval;		
}




