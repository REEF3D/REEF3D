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

#include"wind_f.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice.h"

void wind_f::wind_forcing_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &K, slice &eta)
{
    double windforce,facx,disty;
    
    //cout<<p->mpirank<<"  WIND "<<p->A370<<endl;
    
    if(p->A373<3)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    if(p->A373==1 || eta(i,j)>0.0)
    {
        
    windforce = (p->W3/p->W1)*Cd*p->A371_u*p->A371_u*(cosa*(p->XP[IP]-p->global_xmin) + sina*(p->YP[JP]-p->global_ymin));
    
    K(i,j) -= windforce;
    
    c->test2D(i,j) = windforce;
    
    //cout<<p->mpirank<<" WIND  "<<windforce<<"    "<<p->global_xmin<<endl;
    }
    
    
    if(p->A373==3 || p->A373==4)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    
    {
    Sx = c->Ex(i,j);
    Sy = c->Ey(i,j);
    
    if(Sx*cosa>0.0)
    Sx = 1.0;
    
    else
    Sx = 0.0;
    
    
    if(Sy*sina>0.0)
    Sy = 1.0;
    
    else
    Sy = 0.0;
    
    windforce = 0.0;
    
    if(p->A373==3 || eta(i,j)>0.0)
    windforce = (p->W3/p->W1)*Cd*p->A371_u*p->A371_u*(cosa*(p->XP[IP]-p->global_xmin)*Sx + sina*(p->YP[JP]-p->global_ymin)*Sy);
    
    K(i,j) -= windforce;
    
    c->test2D(i,j) = windforce;
    
    //cout<<p->mpirank<<" WIND  "<<windforce<<"    "<<p->global_xmin<<endl;
    }
    
    LOOP
    c->test(i,j,k) = c->test2D(i,j);
}



